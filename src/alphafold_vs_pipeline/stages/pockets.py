from __future__ import annotations

from typing import Any


def detect_pockets(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    return [{"id": f"P{i+1}", "tool": tool, "score": round(0.9 - i * 0.1, 2)} for i, tool in enumerate(cfg.get("tools", []))]
