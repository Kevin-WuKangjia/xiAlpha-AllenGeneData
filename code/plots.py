from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from config import SingleFeaturePipelineConfig
from utils import zscore_vector


def robust_limits(values: np.ndarray) -> tuple[float, float]:
    """Return adaptive display limits with margin."""
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return -1.0, 1.0
    lo, hi = float(np.min(finite)), float(np.max(finite))
    if np.isclose(lo, hi):
        return lo - 1.0, hi + 1.0
    margin = 0.06 * (hi - lo)
    return lo - margin, hi + margin


def display_mask(x_values: np.ndarray, y_values: np.ndarray) -> np.ndarray:
    """Return an IQR mask that hides extreme points for visualization only."""
    x_values = np.asarray(x_values, dtype=np.float64)
    y_values = np.asarray(y_values, dtype=np.float64)
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if np.sum(finite) < 8:
        return finite

    def fence(values: np.ndarray) -> tuple[float, float]:
        q1, q3 = np.percentile(values, [25.0, 75.0])
        iqr = q3 - q1
        if not np.isfinite(iqr) or iqr <= 0.0:
            return float(np.min(values)), float(np.max(values))
        return float(q1 - 1.5 * iqr), float(q3 + 1.5 * iqr)

    x_lo, x_hi = fence(x_values[finite])
    y_lo, y_hi = fence(y_values[finite])
    mask = finite & (x_values >= x_lo) & (x_values <= x_hi) & (y_values >= y_lo) & (y_values <= y_hi)
    return mask if np.sum(mask) >= 8 else finite


def kde_or_hist(axis: plt.Axes, values: np.ndarray, *, orientation: str, color: str) -> None:
    """Draw a compact marginal density, falling back to a histogram when needed."""
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size < 4 or np.isclose(np.std(values), 0.0):
        if orientation == "x":
            axis.hist(values, bins=12, color=color, alpha=0.70)
        else:
            axis.hist(values, bins=12, orientation="horizontal", color=color, alpha=0.70)
        return
    grid = np.linspace(float(np.min(values)), float(np.max(values)), 160)
    density = gaussian_kde(values)(grid)
    if orientation == "x":
        axis.fill_between(grid, density, color=color, alpha=0.28, linewidth=0)
        axis.plot(grid, density, color=color, linewidth=1.1)
    else:
        axis.fill_betweenx(grid, density, color=color, alpha=0.28, linewidth=0)
        axis.plot(density, grid, color=color, linewidth=1.1)


def plot_pls1_beta2_scatter(config: SingleFeaturePipelineConfig, pls_summary: dict[str, object]) -> Path:
    """Plot z-scored PLS1 score versus z-scored CM beta2 for one feature."""
    score_path = Path(pls_summary["outputs"]["PLS1score"])
    rows = list(csv.DictReader(score_path.open("r", encoding="utf-8-sig", newline="")))
    beta2_raw = np.asarray([float(row["cm_beta2_value"]) for row in rows], dtype=np.float64)
    score_raw = np.asarray([float(row["PLS1_score"]) for row in rows], dtype=np.float64)
    x_values = zscore_vector(beta2_raw)
    y_values = zscore_vector(score_raw)
    mask = display_mask(x_values, y_values)
    x_plot = x_values[mask]
    y_plot = y_values[mask]
    xlim = robust_limits(x_plot)
    ylim = robust_limits(y_plot)

    fig = plt.figure(figsize=(5.4, 5.4), constrained_layout=False)
    ax_main = fig.add_axes([0.18, 0.16, 0.58, 0.58])
    ax_top = fig.add_axes([0.18, 0.76, 0.58, 0.12], sharex=ax_main)
    ax_right = fig.add_axes([0.78, 0.16, 0.14, 0.58], sharey=ax_main)

    ax_main.scatter(
        x_plot,
        y_plot,
        s=42,
        color="#2f5d8c",
        alpha=0.88,
        edgecolor="white",
        linewidth=0.55,
    )
    if x_plot.size >= 2 and np.nanstd(x_plot) > 0:
        slope, intercept = np.polyfit(x_plot, y_plot, deg=1)
        x_line = np.linspace(xlim[0], xlim[1], 100)
        ax_main.plot(x_line, slope * x_line + intercept, color="#b2182b", linewidth=1.7)

    ax_main.text(
        0.04,
        0.96,
        f"r = {float(pls_summary['brain_score_cm_beta2_corr']):.3f}\n"
        f"spin p = {float(pls_summary['spin_p_fixed_pls1_score_corr']):.4f}",
        transform=ax_main.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#d0d0d0", "alpha": 0.9},
    )
    ax_main.set_xlabel("Schaefer100 CM beta2 (z)", fontsize=9)
    ax_main.set_ylabel("PLS1 score (z)", fontsize=9)
    ax_main.set_xlim(*xlim)
    ax_main.set_ylim(*ylim)
    ax_main.grid(True, color="#e8e8e8", linewidth=0.6)
    ax_main.spines[["top", "right"]].set_visible(False)
    ax_main.set_box_aspect(1.0)

    ax_top.set_xlim(*xlim)
    kde_or_hist(ax_top, x_plot, orientation="x", color="#d95f02")
    ax_top.set_axis_off()
    ax_right.set_ylim(*ylim)
    kde_or_hist(ax_right, y_plot, orientation="y", color="#377eb8")
    ax_right.set_axis_off()

    output_path = (
        config.output_root
        / "pls"
        / config.feature
        / "figures"
        / f"{config.feature}_PLS1_score_vs_cm_beta2_zscore.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.suptitle(f"{config.short_name}: PLS1 Score vs CM beta2", fontsize=12, y=0.98)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_path

