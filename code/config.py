from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
WORK_ROOT = PIPELINE_ROOT.parent
DEFAULT_EXAMPLE_ROOT = PIPELINE_ROOT / "example"
DEFAULT_EXAMPLE_DATA_ROOT = DEFAULT_EXAMPLE_ROOT / "data"
DEFAULT_TOOLBOX_ROOT = WORK_ROOT / "XiAlpha_AHBA_toolbox"
DEFAULT_VERTEX_MATRIX_ROOT = DEFAULT_EXAMPLE_DATA_ROOT / "vertex_matrices_processed"
DEFAULT_AHBA_MAT = DEFAULT_EXAMPLE_DATA_ROOT / "AHBA" / "ROIxGene_Schaefer100_INT_zscore.mat"
DEFAULT_PARCELLATION_DIR = DEFAULT_EXAMPLE_DATA_ROOT / "parcellation"
DEFAULT_SPHERE_DIR = DEFAULT_EXAMPLE_DATA_ROOT / "fsaverage"
DEFAULT_OUTPUT_ROOT = DEFAULT_EXAMPLE_ROOT

MAP_TOOLKIT_ROOT_TEXT = os.environ.get("MATLAB_TO_PY_MAPS_ROOT", "")
DEFAULT_MAP_TOOLKIT_ROOT = Path(MAP_TOOLKIT_ROOT_TEXT) if MAP_TOOLKIT_ROOT_TEXT else None


@dataclass(frozen=True)
class SingleFeaturePipelineConfig:
    """Configuration for one beta2-AHBA PLS analysis."""

    feature: str
    short_name: str
    output_root: Path = DEFAULT_OUTPUT_ROOT
    toolbox_root: Path = DEFAULT_TOOLBOX_ROOT
    vertex_matrix_root: Path = DEFAULT_VERTEX_MATRIX_ROOT
    ahba_mat: Path = DEFAULT_AHBA_MAT
    parcellation_dir: Path = DEFAULT_PARCELLATION_DIR
    sphere_dir: Path | None = DEFAULT_SPHERE_DIR
    map_toolkit_root: Path | None = DEFAULT_MAP_TOOLKIT_ROOT
    parcel_count: int = 100
    n_bootstrap: int = 1000
    n_spins: int = 5000
    random_state: int = 20260706
    run_zig: bool = True
    existing_zig_root: Path | None = None
