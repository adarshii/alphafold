from __future__ import annotations

from typing import Any


def prepare_structure(target_cfg: dict[str, Any]) -> dict[str, Any]:
    return {
        "source": target_cfg.get("structure_source", "provided"),
        "predictor": target_cfg.get("predictor", "colabfold"),
        "structure_path": target_cfg.get("structure_path", "data/raw/6LU7.pdb"),
        "supported_predictors": ["alphafold2", "esmfold", "colabfold"],
    }
