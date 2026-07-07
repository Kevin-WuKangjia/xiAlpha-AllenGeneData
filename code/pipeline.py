from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from config import SingleFeaturePipelineConfig
from phenotype import parcellate_cm_beta2
from plots import plot_pls1_beta2_scatter
from pls import run_pls_feature
from utils import write_json
from zig_cm import run_zig_cm_feature, zig_cm_figure_path


def run_single_feature_pipeline(config: SingleFeaturePipelineConfig) -> dict[str, object]:
    """Run the full single-feature ZIG/CM to AHBA PLS pipeline."""
    config.output_root.mkdir(parents=True, exist_ok=True)
    zig_root = config.existing_zig_root
    zig_summary = None
    if config.run_zig:
        zig_summary = run_zig_cm_feature(config)
        zig_root = config.output_root / "zig_cm"
    if zig_root is None:
        raise ValueError("existing_zig_root must be provided when run_zig is False.")

    phenotype = parcellate_cm_beta2(config, zig_root=zig_root)
    pls_summary = run_pls_feature(config, phenotype)
    scatter = plot_pls1_beta2_scatter(config, pls_summary)
    zig_figure = zig_cm_figure_path(config, zig_root)

    pipeline_summary = {
        "config": {key: str(value) if isinstance(value, Path) else value for key, value in asdict(config).items()},
        "zig_cm": zig_summary if zig_summary is not None else {"used_existing_zig_root": str(zig_root)},
        "phenotype": phenotype["summary"],
        "pls": pls_summary,
        "figures": {
            "zig_cm_direction_maps": str(zig_figure) if zig_figure.exists() else None,
            "pls1_beta2_scatter_zscore": str(scatter),
        },
        "output_root": str(config.output_root),
    }
    write_json(config.output_root / "summary" / "pipeline_summary.json", pipeline_summary)
    return pipeline_summary

