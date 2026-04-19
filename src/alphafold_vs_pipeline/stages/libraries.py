from __future__ import annotations

from typing import Any


def prepare_library(_cfg: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {"compound_id": "CHEMBL_A1", "smiles": "CCO", "mw": 46.1, "logp": -0.3, "hbd": 1, "hba": 1, "pains": False},
        {"compound_id": "CHEMBL_A2", "smiles": "CCN(CC)CC", "mw": 101.2, "logp": 1.3, "hbd": 0, "hba": 1, "pains": False},
        {"compound_id": "ZINC_D1", "smiles": "c1ccccc1", "mw": 78.1, "logp": 2.1, "hbd": 0, "hba": 0, "pains": False},
    ]
