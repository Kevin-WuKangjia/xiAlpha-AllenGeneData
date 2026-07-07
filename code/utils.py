from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
TOOLBOX_ROOT = PIPELINE_ROOT.parent / "XiAlpha_AHBA_toolbox"
MAP_TOOLKIT_ROOT_TEXT = os.environ.get("MATLAB_TO_PY_MAPS_ROOT", "")
MAP_TOOLKIT_ROOT = Path(MAP_TOOLKIT_ROOT_TEXT) if MAP_TOOLKIT_ROOT_TEXT else None


def configure_import_paths(toolbox_root: Path, map_toolkit_root: Path | None = None) -> None:
    """Add toolbox code folders to sys.path."""
    paths = [
        toolbox_root / "code",
        toolbox_root / "test" / "code",
    ]
    toolkit_text = os.environ.get("MATLAB_TO_PY_MAPS_ROOT", "")
    if map_toolkit_root is not None:
        paths.append(map_toolkit_root)
    elif toolkit_text:
        paths.append(Path(toolkit_text))
    for path in paths:
        if path.exists() and str(path) not in sys.path:
            sys.path.insert(0, str(path))


def compact_local_paths(payload: object) -> object:
    """Replace machine-specific local roots in JSON payloads with portable labels."""
    replacements = [
        (str(PIPELINE_ROOT), "<PIPELINE_ROOT>"),
        (str(TOOLBOX_ROOT), "<TOOLBOX_ROOT>"),
    ]
    if MAP_TOOLKIT_ROOT is not None:
        replacements.append((str(MAP_TOOLKIT_ROOT), "<MAP_TOOLKIT_ROOT>"))

    def compact_value(value: object) -> object:
        if isinstance(value, dict):
            return {key: compact_value(item) for key, item in value.items()}
        if isinstance(value, list):
            return [compact_value(item) for item in value]
        if isinstance(value, tuple):
            return [compact_value(item) for item in value]
        if isinstance(value, str):
            compacted = value
            for root, label in replacements:
                compacted = compacted.replace(root, label)
                compacted = compacted.replace(root.replace("\\", "/"), label)
            return compacted.replace("\\", "/")
        return value

    return compact_value(payload)


def write_json(path: Path, payload: object) -> None:
    """Write one JSON file with UTF-8 encoding."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(compact_local_paths(payload), indent=2, ensure_ascii=False), encoding="utf-8")


def write_matrix_csv(path: Path, row_name: str, column_names: list[str], values: np.ndarray) -> None:
    """Write a one-row labeled matrix CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    values = np.asarray(values, dtype=np.float64).reshape(1, -1)
    with path.open("w", encoding="utf-8", newline="") as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(["parameter", *column_names])
        writer.writerow([row_name, *[f"{value:.12g}" for value in values[0]]])


def zscore_vector(values: np.ndarray) -> np.ndarray:
    """Z-score one vector, returning zeros if the vector has no variance."""
    values = np.asarray(values, dtype=np.float64)
    std_value = float(np.nanstd(values))
    if not np.isfinite(std_value) or std_value <= 0.0:
        return np.zeros_like(values, dtype=np.float64)
    return (values - float(np.nanmean(values))) / std_value
