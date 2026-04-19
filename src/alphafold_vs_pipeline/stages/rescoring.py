from __future__ import annotations

from typing import Any


def rescore_poses(compounds: list[dict[str, Any]], _cfg: dict[str, Any]) -> list[dict[str, Any]]:
    for item in compounds:
        item["rf_score_features"] = {"hbonds": 1, "hydrophobics": 2, "metal_contacts": 0}
        item["ml_rescore"] = round(item["docking_score"] * 1.15, 3)
    return compounds
