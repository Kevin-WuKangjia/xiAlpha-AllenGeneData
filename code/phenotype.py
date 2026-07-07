from __future__ import annotations

from pathlib import Path

import numpy as np

from config import SingleFeaturePipelineConfig
from utils import configure_import_paths, write_json, write_matrix_csv
from zig_cm import cm_beta2_paths


def parcellate_cm_beta2(config: SingleFeaturePipelineConfig, zig_root: Path | None = None) -> dict[str, object]:
    """Parcellate bilateral CM beta2 maps to Schaefer100."""
    configure_import_paths(config.toolbox_root, config.map_toolkit_root)
    from parcellate_fsaverage10k import parcellate_fsaverage10k

    left_gii, right_gii = cm_beta2_paths(config, zig_root)
    region_names, bilateral_values = parcellate_fsaverage10k(
        [left_gii, right_gii],
        hemi="BOTH",
        parcel_count=config.parcel_count,
        parcellation_dir=config.parcellation_dir,
        sphere_dir=config.sphere_dir,
    )
    bilateral_values = np.asarray(bilateral_values, dtype=np.float64)
    left_count = config.parcel_count // 2
    left_names = region_names[:left_count]
    left_values = bilateral_values[:left_count]

    output_dir = config.output_root / "phenotype"
    write_matrix_csv(output_dir / "cm_beta2_bilateral.csv", config.feature, region_names, bilateral_values)
    write_matrix_csv(output_dir / "cm_beta2_left.csv", config.feature, left_names, left_values)
    np.save(output_dir / "cm_beta2_bilateral.npy", bilateral_values.reshape(1, -1))
    np.save(output_dir / "cm_beta2_left.npy", left_values.reshape(1, -1))

    summary = {
        "module": "phenotype",
        "feature": config.feature,
        "input_left_gii": str(left_gii),
        "input_right_gii": str(right_gii),
        "parcel_count": config.parcel_count,
        "left_roi_count": left_count,
        "outputs": {
            "bilateral_csv": str(output_dir / "cm_beta2_bilateral.csv"),
            "bilateral_npy": str(output_dir / "cm_beta2_bilateral.npy"),
            "left_csv": str(output_dir / "cm_beta2_left.csv"),
            "left_npy": str(output_dir / "cm_beta2_left.npy"),
        },
    }
    write_json(config.output_root / "summary" / "phenotype_summary.json", summary)
    return {
        "summary": summary,
        "region_names_bilateral": region_names,
        "region_names_left": left_names,
        "values_bilateral": bilateral_values,
        "values_left": left_values,
    }
