from __future__ import annotations

import hashlib
from typing import Any


def _deterministic_score(compound_id: str, base: float = -6.0, spread: float = 4.0) -> float:
    digest = hashlib.sha256(compound_id.encode("utf-8")).hexdigest()
    fraction = int(digest[:8], 16) / 0xFFFFFFFF
    return round(base - (fraction * spread), 3)


def dock_batch(compounds: list[dict[str, Any]], _cfg: dict[str, Any]) -> list[dict[str, Any]]:
    for item in compounds:
        item["docking_score"] = _deterministic_score(item["compound_id"])
    return compounds
