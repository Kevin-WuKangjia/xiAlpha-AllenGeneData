from __future__ import annotations

import argparse
import sys
from pathlib import Path


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = PIPELINE_ROOT / "code"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from config import SingleFeaturePipelineConfig
from pipeline import run_single_feature_pipeline


DEFAULT_OUTPUT_ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    """Parse command-line options for the Alpha Power example."""
    parser = argparse.ArgumentParser(description="Run the Alpha Power beta2-AHBA PLS example.")
    parser.add_argument(
        "--skip-zig",
        action="store_true",
        help="Use an existing ZIG/CM output directory instead of refitting ZIG/CM.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Example output directory.",
    )
    parser.add_argument(
        "--existing-zig-root",
        type=Path,
        default=None,
        help="Existing ZIG/CM root used with --skip-zig. Defaults to <output-root>/zig_cm.",
    )
    parser.add_argument("--n-bootstrap", type=int, default=1000, help="Number of ROI bootstraps for gene stability.")
    parser.add_argument("--n-spins", type=int, default=5000, help="Number of spin permutations.")
    return parser.parse_args()


def main() -> None:
    """Run the single-feature pipeline on Alpha_estimate_Power."""
    args = parse_args()
    existing_zig_root = None
    if args.skip_zig:
        existing_zig_root = args.existing_zig_root or args.output_root / "zig_cm"

    config = SingleFeaturePipelineConfig(
        feature="Alpha_estimate_Power",
        short_name="Aalpha",
        output_root=args.output_root,
        n_bootstrap=args.n_bootstrap,
        n_spins=args.n_spins,
        random_state=20260706,
        run_zig=not args.skip_zig,
        existing_zig_root=existing_zig_root,
    )
    summary = run_single_feature_pipeline(config)
    print(f"[DONE] {summary['output_root']}")
    print(f"[SUMMARY] {args.output_root / 'summary' / 'pipeline_summary.json'}")


if __name__ == "__main__":
    main()
