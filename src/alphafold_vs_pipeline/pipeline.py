from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from alphafold_vs_pipeline.stages import (
    detect_pockets,
    dock_batch,
    filter_admet,
    prepare_library,
    prepare_structure,
    rank_and_export,
    rescore_poses,
)


def run_pipeline(config: dict[str, Any], output_dir: str | Path, dry_run: bool = False) -> dict[str, Any]:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    structure_info = prepare_structure(config["target"])
    pockets = detect_pockets(config["pocket_detection"])
    compounds = prepare_library(config["library"])
    docked = dock_batch(compounds, config["docking"])
    rescored = rescore_poses(docked, config["rescoring"])
    passed = filter_admet(rescored, config["admet"])
    ranked = rank_and_export(passed, config["ranking"], out_dir)

    summary = {
        "target": config["target"]["name"],
        "target_pdb_id": config["target"]["pdb_id"],
        "structure": structure_info,
        "pockets": pockets,
        "hits": ranked,
        "mode": "dry-run" if dry_run else "standard",
    }

    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    return summary
