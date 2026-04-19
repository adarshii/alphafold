from __future__ import annotations

import argparse
import json
from pathlib import Path

from alphafold_vs_pipeline.config import load_config
from alphafold_vs_pipeline.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="alphafold-vs", description="AlphaFold-guided virtual screening pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run end-to-end VS pipeline")
    run_parser.add_argument("--config", required=True, help="Path to YAML config file")
    run_parser.add_argument("--output", required=True, help="Output directory")
    run_parser.add_argument("--dry-run", action="store_true", help="Run mock/demo pipeline without external tools")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        config = load_config(args.config)
        result = run_pipeline(config=config, output_dir=Path(args.output), dry_run=args.dry_run)
        print(json.dumps(result, indent=2))
        return 0
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
