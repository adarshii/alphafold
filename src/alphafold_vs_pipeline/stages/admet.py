from __future__ import annotations

from typing import Any


def filter_admet(compounds: list[dict[str, Any]], cfg: dict[str, Any]) -> list[dict[str, Any]]:
    filtered: list[dict[str, Any]] = []
    for item in compounds:
        if item["mw"] > cfg.get("max_mw", 500):
            continue
        if item["logp"] > cfg.get("max_logp", 5):
            continue
        if item["hbd"] > cfg.get("max_hbd", 5):
            continue
        if item["hba"] > cfg.get("max_hba", 10):
            continue
        if cfg.get("pains_filter", True) and item["pains"]:
            continue
        filtered.append(item)
    return filtered
