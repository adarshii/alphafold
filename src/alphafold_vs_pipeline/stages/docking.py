from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


def _deterministic_score(compound_id: str, base: float = -6.0, spread: float = 4.0) -> float:
    digest = hashlib.sha256(compound_id.encode("utf-8")).hexdigest()
    fraction = int(digest[:8], 16) / 0xFFFFFFFF
    return round(base - (fraction * spread), 3)


def _prepare_ligand_pdbqt(smiles: str, ligand_id: str, out_dir: Path) -> Path:
    sdf_path = out_dir / f"{ligand_id}.sdf"
    pdbqt_path = out_dir / f"{ligand_id}.pdbqt"

    sdf_path.write_text(f"{smiles}\n", encoding="utf-8")
    command = ["obabel", "-isdf", str(sdf_path), "-opdbqt", "-O", str(pdbqt_path), "--gen3d"]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return pdbqt_path


def _prepare_protein_pdbqt(structure_path: str, out_dir: Path) -> Path:
    receptor_pdbqt = out_dir / "receptor.pdbqt"
    command = ["obabel", "-ipdb", structure_path, "-opdbqt", "-O", str(receptor_pdbqt)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return receptor_pdbqt


def _run_vina(
    executable: str,
    receptor: Path,
    ligand: Path,
    center: list[float],
    size: list[float],
    out_pose: Path,
    cfg: dict[str, Any],
) -> float:
    command = [
        executable,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--exhaustiveness",
        str(cfg.get("exhaustiveness", 8)),
        "--num_modes",
        str(cfg.get("num_modes", 10)),
        "--cpu",
        str(cfg.get("cpu", 1)),
        "--out",
        str(out_pose),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("1 "):
            parts = stripped.split()
            return float(parts[1])
    raise RuntimeError("Unable to parse Vina affinity from output")


def dock_batch(
    compounds: list[dict[str, Any]],
    cfg: dict[str, Any],
    structure_info: dict[str, Any],
    pockets: list[dict[str, Any]],
    out_dir: Path,
    dry_run: bool = False,
    seed: int = 42,
) -> list[dict[str, Any]]:
    docking_dir = out_dir / "docking"
    docking_dir.mkdir(parents=True, exist_ok=True)

    if not pockets:
        raise ValueError("No pocket metadata available for docking")
    pocket = pockets[0]

    if dry_run:
        for item in compounds:
            item["docking_score"] = _deterministic_score(item["compound_id"])
            item["pose_path"] = str(docking_dir / f"{item['compound_id']}_pose.pdbqt")
        (docking_dir / "docking_results.json").write_text(json.dumps(compounds, indent=2), encoding="utf-8")
        return compounds

    executable = str(cfg.get("executable", cfg.get("engine", "vina")))
    protein_pdbqt = cfg.get("protein_pdbqt")
    receptor = Path(protein_pdbqt) if protein_pdbqt else _prepare_protein_pdbqt(structure_info["prepared_structure_path"], docking_dir)

    for item in compounds:
        ligand_id = str(item["compound_id"])
        ligand_pdbqt = _prepare_ligand_pdbqt(item["smiles"], ligand_id, docking_dir)
        pose_path = docking_dir / f"{ligand_id}_pose.pdbqt"
        affinity = _run_vina(
            executable=executable,
            receptor=receptor,
            ligand=ligand_pdbqt,
            center=pocket["center"],
            size=pocket["size"],
            out_pose=pose_path,
            cfg=cfg,
        )
        item["docking_score"] = affinity
        item["pose_path"] = str(pose_path)

    (docking_dir / "docking_results.json").write_text(json.dumps(compounds, indent=2), encoding="utf-8")
    return compounds
