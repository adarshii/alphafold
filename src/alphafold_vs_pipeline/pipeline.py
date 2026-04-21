from __future__ import annotations

import json
import logging
import random
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from alphafold_vs_pipeline.stages import (
    apply_admet,
    detect_pockets,
    dock_batch,
    prepare_library,
    prepare_structure,
    rank_and_export,
    rescore_poses,
)

LOGGER = logging.getLogger(__name__)


def _set_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def run_pipeline(
    config: dict[str, Any],
    output_dir: str | Path,
    dry_run: bool = False,
    stage: str = "full",
) -> dict[str, Any]:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    reproducibility_cfg = config.get("reproducibility", {})
    seed = int(reproducibility_cfg.get("random_seed", 42))
    _set_global_seed(seed)

    (out_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    structure_info: dict[str, Any] = {}
    pockets: list[dict[str, Any]] = []
    compounds: list[dict[str, Any]] = []

    structure_info = prepare_structure(config["target"], out_dir=out_dir, dry_run=dry_run)
    if stage == "structure":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed)

    pockets = detect_pockets(config["pocket_detection"], structure_info=structure_info, out_dir=out_dir, dry_run=dry_run)
    if stage == "pockets":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed)

    compounds = prepare_library(config["library"], out_dir=out_dir, dry_run=dry_run, seed=seed)
    if stage == "library":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed)

    compounds = dock_batch(
        compounds,
        config["docking"],
        structure_info=structure_info,
        pockets=pockets,
        out_dir=out_dir,
        dry_run=dry_run,
        seed=seed,
    )
    if stage == "docking":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed)

    compounds, ml_artifacts = rescore_poses(
        compounds,
        config["rescoring"],
        benchmark_cfg=config.get("benchmark", {}),
        out_dir=out_dir,
        dry_run=dry_run,
        seed=seed,
    )
    if stage == "rescoring":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed, ml_artifacts)

    compounds = apply_admet(compounds, config["admet"], out_dir=out_dir)
    if stage == "admet":
        return _write_summary(config, out_dir, structure_info, pockets, compounds, dry_run, seed, ml_artifacts)

    ranked = rank_and_export(compounds, config["ranking"], out_dir)
    return _write_summary(config, out_dir, structure_info, pockets, ranked, dry_run, seed, ml_artifacts)


def _write_summary(
    config: dict[str, Any],
    out_dir: Path,
    structure_info: dict[str, Any],
    pockets: list[dict[str, Any]],
    compounds: list[dict[str, Any]],
    dry_run: bool,
    seed: int,
    ml_artifacts: dict[str, Any] | None = None,
) -> dict[str, Any]:
    summary = {
        "target": config["target"]["name"],
        "target_pdb_id": config["target"].get("pdb_id", "N/A"),
        "structure": structure_info,
        "pockets": pockets,
        "hits": compounds,
        "mode": "dry-run" if dry_run else "standard",
        "seed": seed,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "ml": ml_artifacts or {},
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    LOGGER.info("Summary Written to %s", summary_path)
    return summary
