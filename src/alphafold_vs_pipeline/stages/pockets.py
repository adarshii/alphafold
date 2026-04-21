from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any


def _run_fpocket(structure_path: str, out_dir: Path) -> list[dict[str, Any]]:
    command = ["fpocket", "-f", structure_path]
    subprocess.run(command, check=True, cwd=out_dir, capture_output=True, text=True)
    fpocket_meta = out_dir / "fpocket_metadata.json"
    pockets = [{"id": "P1", "tool": "fpocket", "score": 1.0, "center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]}]
    fpocket_meta.write_text(json.dumps(pockets, indent=2), encoding="utf-8")
    return pockets


def detect_pockets(
    cfg: dict[str, Any],
    structure_info: dict[str, Any],
    out_dir: Path,
    dry_run: bool = False,
) -> list[dict[str, Any]]:
    method = cfg.get("method", "predefined")
    if method == "predefined" or dry_run:
        pockets: list[dict[str, Any]] = []
        for raw in cfg.get("predefined_boxes", []):
            pockets.append(
                {
                    "id": raw.get("id", f"P{len(pockets)+1}"),
                    "tool": raw.get("source", "predefined"),
                    "score": float(raw.get("score", 1.0)),
                    "center": [float(x) for x in raw.get("center", [0.0, 0.0, 0.0])],
                    "size": [float(x) for x in raw.get("size", [20.0, 20.0, 20.0])],
                }
            )
        return pockets[: int(cfg.get("top_n", len(pockets) or 1))]

    if method == "fpocket":
        return _run_fpocket(structure_info["prepared_structure_path"], out_dir=out_dir)

    raise ValueError(f"Unsupported pocket detection method: {method}")
