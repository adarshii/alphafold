from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any


def prepare_structure(target_cfg: dict[str, Any], out_dir: Path, dry_run: bool = False) -> dict[str, Any]:
    structure_path = Path(target_cfg.get("structure_path", ""))
    if not dry_run and not structure_path.exists():
        raise FileNotFoundError(f"Structure path does not exist: {structure_path}")

    copied_path = out_dir / "target_structure.pdb"
    if not dry_run and structure_path.exists():
        shutil.copyfile(structure_path, copied_path)

    return {
        "source": target_cfg.get("structure_source", "provided"),
        "predictor": target_cfg.get("predictor", "colabfold"),
        "structure_path": str(structure_path),
        "prepared_structure_path": str(copied_path),
        "supported_predictors": ["alphafold2", "esmfold", "colabfold"],
    }
