from __future__ import annotations

from pathlib import Path

from config import SingleFeaturePipelineConfig
from utils import configure_import_paths, write_json


def run_zig_cm_feature(config: SingleFeaturePipelineConfig) -> dict[str, object]:
    """Run vertex-wise ZIP gamma1 and CM beta2 models for one feature."""
    configure_import_paths(config.toolbox_root, config.map_toolkit_root)
    import run_vertex_zig_zip_cm_all as zig

    zig.MATRIX_ROOT = config.vertex_matrix_root
    zig.OUTPUT_ROOT = config.output_root / "zig_cm"
    summary = zig.process_feature(config.feature, config.short_name)
    module_summary = {
        "module": "zig_cm",
        "feature": config.feature,
        "short_name": config.short_name,
        "input_vertex_matrix_root": str(config.vertex_matrix_root),
        "output_root": str(zig.OUTPUT_ROOT),
        "model": {
            "zip": "logit(P(Y == 0)) = gamma0 + gamma1 * z_age",
            "cm": "Y_nonzero = beta0 + beta1 * z_age + beta2 * z_age^2",
        },
        "feature_summary": summary,
    }
    write_json(config.output_root / "summary" / "zig_cm_summary.json", module_summary)
    return module_summary


def cm_beta2_paths(config: SingleFeaturePipelineConfig, zig_root: Path | None = None) -> tuple[Path, Path]:
    """Return left and right raw CM beta2 GIFTI paths."""
    root = zig_root or config.output_root / "zig_cm"
    feature_dir = root / config.feature
    left = feature_dir / f"{config.feature}_cm_beta2_raw_hemi-L.shape.gii"
    right = feature_dir / f"{config.feature}_cm_beta2_raw_hemi-R.shape.gii"
    if not left.exists() or not right.exists():
        raise FileNotFoundError(f"Missing CM beta2 GIFTI files: {left}, {right}")
    return left, right


def zig_cm_figure_path(config: SingleFeaturePipelineConfig, zig_root: Path) -> Path:
    """Return the expected ZIG/CM direction-map figure path."""
    return zig_root / config.feature / f"{config.feature}_ZIG_ZIP_CM_direction_2x4.png"
