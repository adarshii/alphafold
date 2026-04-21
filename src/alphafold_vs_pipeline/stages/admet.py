from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _toxicity_flags(smiles: str) -> dict[str, bool]:
    try:
        from rdkit import Chem
    except ImportError:
        lower = smiles.lower()
        return {
            "nitro_alert": "[n+](=o)[o-]" in lower,
            "aniline_alert": "nc1ccccc1" in lower,
            "highly_halogenated": sum(lower.count(x) for x in ["cl", "br", "i"]) >= 4,
        }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"nitro_alert": False, "aniline_alert": False, "highly_halogenated": False}

    nitro = Chem.MolFromSmarts("[NX3+](=O)[O-]")
    aniline = Chem.MolFromSmarts("Nc1ccccc1")
    halogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in {"Cl", "Br", "I"})
    return {
        "nitro_alert": bool(nitro and mol.HasSubstructMatch(nitro)),
        "aniline_alert": bool(aniline and mol.HasSubstructMatch(aniline)),
        "highly_halogenated": halogen_atoms >= 4,
    }


def _drug_likeness_score(item: dict[str, Any], thresholds: dict[str, float], toxicity: dict[str, bool]) -> float:
    checks = [
        item["mw"] <= thresholds["max_mw"],
        item["logp"] <= thresholds["max_logp"],
        item["hbd"] <= thresholds["max_hbd"],
        item["hba"] <= thresholds["max_hba"],
        item["tpsa"] <= thresholds["max_tpsa"],
        item["rotatable_bonds"] <= thresholds["max_rotatable_bonds"],
    ]
    toxicity_penalty = sum(1 for v in toxicity.values() if v)
    return float(max(0.0, (sum(checks) / len(checks)) - 0.1 * toxicity_penalty))


def apply_admet(compounds: list[dict[str, Any]], cfg: dict[str, Any], out_dir: Path) -> list[dict[str, Any]]:
    thresholds = {
        "max_mw": float(cfg.get("max_mw", 500)),
        "max_logp": float(cfg.get("max_logp", 5)),
        "max_hbd": float(cfg.get("max_hbd", 5)),
        "max_hba": float(cfg.get("max_hba", 10)),
        "max_tpsa": float(cfg.get("max_tpsa", 140)),
        "max_rotatable_bonds": float(cfg.get("max_rotatable_bonds", 10)),
    }

    max_violations = int(cfg.get("max_violations", 1))
    apply_filters = bool(cfg.get("apply_filters", True))

    profiled: list[dict[str, Any]] = []
    for item in compounds:
        toxicity = _toxicity_flags(str(item["smiles"]))
        violations = sum(
            [
                int(item["mw"] > thresholds["max_mw"]),
                int(item["logp"] > thresholds["max_logp"]),
                int(item["hbd"] > thresholds["max_hbd"]),
                int(item["hba"] > thresholds["max_hba"]),
                int(item["tpsa"] > thresholds["max_tpsa"]),
                int(item["rotatable_bonds"] > thresholds["max_rotatable_bonds"]),
            ]
        )

        item["toxicity_flags"] = toxicity
        item["admet_violations"] = violations
        item["drug_likeness_score"] = _drug_likeness_score(item, thresholds, toxicity)

        if apply_filters and violations > max_violations:
            continue
        profiled.append(item)

    (out_dir / "admet_profile.json").write_text(json.dumps(profiled, indent=2), encoding="utf-8")
    return profiled
