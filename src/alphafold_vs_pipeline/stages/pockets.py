from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any


def _parse_atom_coordinates(pdb_path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            coords.append((x, y, z))
    return coords


def _center_and_size(coords: list[tuple[float, float, float]]) -> tuple[list[float], list[float]]:
    xs = [p[0] for p in coords]
    ys = [p[1] for p in coords]
    zs = [p[2] for p in coords]
    center = [float((min(xs) + max(xs)) / 2), float((min(ys) + max(ys)) / 2), float((min(zs) + max(zs)) / 2)]
    size = [float(max(8.0, max(xs) - min(xs) + 4.0)), float(max(8.0, max(ys) - min(ys) + 4.0)), float(max(8.0, max(zs) - min(zs) + 4.0))]
    return center, size


def _run_fpocket(structure_path: str, out_dir: Path) -> list[dict[str, Any]]:
    command = ["fpocket", "-f", structure_path]
    try:
        subprocess.run(command, check=True, cwd=out_dir, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError("fpocket command failed: fpocket is not installed or not in PATH.") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"fpocket command failed: {exc.stderr}\n{exc.stdout}") from exc

    structure_stem = Path(structure_path).stem
    fpocket_out = out_dir / f"{structure_stem}_out" / "pockets"
    pocket_files = sorted(fpocket_out.glob("pocket*_atm.pdb"))
    if not pocket_files:
        raise RuntimeError(f"fpocket completed but no pocket files were found in {fpocket_out}")

    fpocket_meta = out_dir / "fpocket_metadata.json"
    pockets: list[dict[str, Any]] = []
    for index, pocket_path in enumerate(pocket_files, start=1):
        coords = _parse_atom_coordinates(pocket_path)
        if not coords:
            continue
        center, size = _center_and_size(coords)
        pockets.append(
            {
                "id": f"P{index}",
                "tool": "fpocket",
                "score": float(len(coords)),
                "center": center,
                "size": size,
                "source_file": str(pocket_path),
            }
        )
    if not pockets:
        raise RuntimeError("fpocket output was found but no valid pocket coordinates could be parsed.")
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
