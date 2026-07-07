from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from config import SingleFeaturePipelineConfig
from utils import configure_import_paths, write_json, zscore_vector


def split_genes_by_bootstrap_ratio(gene_weights_csv: Path, output_dir: Path, threshold: float = 3.0) -> dict[str, object]:
    """Save positive and negative gene lists by bootstrap-ratio threshold in descending order."""
    with gene_weights_csv.open("r", encoding="utf-8-sig", newline="") as f_obj:
        rows = list(csv.DictReader(f_obj))
    positive = [row for row in rows if float(row["bootstrap_ratio"]) > threshold]
    negative = [row for row in rows if float(row["bootstrap_ratio"]) < -threshold]
    positive.sort(key=lambda row: float(row["bootstrap_ratio"]), reverse=True)
    negative.sort(key=lambda row: float(row["bootstrap_ratio"]), reverse=True)
    outputs: dict[str, object] = {
        "threshold": threshold,
        "positive_gene_count": len(positive),
        "negative_gene_count": len(negative),
    }
    for label, selected in (("positive", positive), ("negative", negative)):
        path = output_dir / f"genes_bootstrap_ratio_{label}_gt{threshold:g}.csv"
        with path.open("w", encoding="utf-8", newline="") as f_obj:
            writer = csv.writer(f_obj)
            writer.writerow(["gene", "bootstrap_ratio"])
            for row in selected:
                writer.writerow([row["gene"], row["bootstrap_ratio"]])
        outputs[f"{label}_genes"] = str(path)
    return outputs


def run_fixed_pls1_spin_test(
    *,
    brain_scores: np.ndarray,
    bilateral_values: np.ndarray,
    spin_indices: np.ndarray,
    output_dir: Path,
) -> dict[str, object]:
    """Spin the bilateral phenotype and correlate spun left maps with the fixed real PLS1 score."""
    left_count = brain_scores.size
    real_corr = float(np.corrcoef(brain_scores, bilateral_values[:left_count])[0, 1])
    null_corr = np.zeros(spin_indices.shape[0], dtype=np.float64)
    for spin_idx, spin_index in enumerate(spin_indices):
        spun_left = bilateral_values[spin_index][:left_count]
        null_corr[spin_idx] = float(np.corrcoef(brain_scores, spun_left)[0, 1])
    p_value = float((1 + np.sum(null_corr >= real_corr)) / (null_corr.size + 1))

    null_path = output_dir / "spin_null_fixed_pls1_score_corr.csv"
    with null_path.open("w", encoding="utf-8", newline="") as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(["spin_index", "corr"])
        for idx, value in enumerate(null_corr):
            writer.writerow([idx, f"{value:.12g}"])

    summary = {
        "n_spins": int(null_corr.size),
        "real_pls1_score_cm_beta2_corr": real_corr,
        "spin_p_fixed_pls1_score_corr": p_value,
        "null_mean": float(np.mean(null_corr)),
        "null_std": float(np.std(null_corr)),
        "null_max": float(np.max(null_corr)),
        "test_tail": "one-sided greater-or-equal",
        "null_distribution_csv": str(null_path),
    }
    write_json(output_dir / "spin_test_fixed_pls1_summary.json", summary)
    return summary


def write_pls1score_table(
    output_dir: Path,
    region_names: list[str],
    cm_beta2_values: np.ndarray,
    brain_scores: np.ndarray,
) -> Path:
    """Write the ROI-wise phenotype and PLS1 score table used for plotting."""
    path = output_dir / "PLS1score.csv"
    beta2_z = zscore_vector(cm_beta2_values)
    with path.open("w", encoding="utf-8", newline="") as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(["region", "cm_beta2_value", "cm_beta2_z", "PLS1_score"])
        for region, beta2, beta2_std, score in zip(region_names, cm_beta2_values, beta2_z, brain_scores, strict=True):
            writer.writerow([region, f"{beta2:.12g}", f"{beta2_std:.12g}", f"{score:.12g}"])
    return path


def run_pls_feature(config: SingleFeaturePipelineConfig, phenotype: dict[str, object]) -> dict[str, object]:
    """Run AHBA PLS and fixed-score spin test for one CM beta2 phenotype."""
    configure_import_paths(config.toolbox_root, config.map_toolkit_root)
    from pls_gene_phenotype import fit_gene_phenotype_pls, load_ahba_gene_matrix, write_pls_result
    from spin_test_pls import make_bilateral_spin_indices, parcel_sphere_centroids

    gene_expression, gene_names, _ = load_ahba_gene_matrix(config.ahba_mat)
    left_values = np.asarray(phenotype["values_left"], dtype=np.float64)
    bilateral_values = np.asarray(phenotype["values_bilateral"], dtype=np.float64)
    region_names_left = list(phenotype["region_names_left"])
    if gene_expression.shape[0] != left_values.size:
        raise ValueError(
            "AHBA rows must match left-hemisphere phenotype ROIs: "
            f"AHBA={gene_expression.shape[0]}, phenotype={left_values.size}"
        )

    output_dir = config.output_root / "pls" / config.feature
    output_dir.mkdir(parents=True, exist_ok=True)
    result = fit_gene_phenotype_pls(
        gene_expression,
        left_values,
        gene_names=gene_names,
        region_names=region_names_left,
        n_bootstrap=config.n_bootstrap,
        random_state=config.random_state,
        standardize=True,
    )
    outputs = write_pls_result(result, output_dir)
    pls1score = write_pls1score_table(output_dir, region_names_left, left_values, result.brain_scores)
    split_outputs = split_genes_by_bootstrap_ratio(Path(outputs["gene_weights"]), output_dir)

    centroids = parcel_sphere_centroids(
        config.parcel_count,
        parcellation_dir=config.parcellation_dir,
        sphere_dir=config.sphere_dir,
    )
    spin_indices = make_bilateral_spin_indices(centroids, n_spins=config.n_spins, random_state=config.random_state)
    spin_summary = run_fixed_pls1_spin_test(
        brain_scores=result.brain_scores,
        bilateral_values=bilateral_values,
        spin_indices=spin_indices,
        output_dir=output_dir,
    )

    summary = {
        "module": "pls",
        "feature": config.feature,
        "short_name": config.short_name,
        "ahba_mat": str(config.ahba_mat),
        "n_bootstrap": config.n_bootstrap,
        "n_spins": config.n_spins,
        "brain_score_cm_beta2_corr": spin_summary["real_pls1_score_cm_beta2_corr"],
        "spin_p_fixed_pls1_score_corr": spin_summary["spin_p_fixed_pls1_score_corr"],
        "singular_value": result.summary["singular_value"],
        "positive_br_gt3_gene_count": split_outputs["positive_gene_count"],
        "negative_br_gt3_gene_count": split_outputs["negative_gene_count"],
        "outputs": {
            **{f"pls_{key}": value for key, value in outputs.items()},
            "PLS1score": str(pls1score),
            **split_outputs,
            "spin_summary_json": str(output_dir / "spin_test_fixed_pls1_summary.json"),
            "spin_null_fixed_pls1_score_corr": str(output_dir / "spin_null_fixed_pls1_score_corr.csv"),
        },
    }
    write_json(config.output_root / "summary" / "pls_summary.json", summary)
    return summary
