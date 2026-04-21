from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


def _safe_import_rdkit() -> Any:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required for ligand preparation; install with full extras.") from exc
    return Chem, AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors


def _fingerprint_bits(all_chem: Any, mol: Any, radius: int = 2, n_bits: int = 1024) -> list[int]:
    fp = all_chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return [int(bit) for bit in fp.ToBitString()]


def _compute_descriptors(descriptors: Any, crippen: Any, lipinski: Any, rd_mol_descriptors: Any, mol: Any) -> dict[str, float]:
    return {
        "mw": float(descriptors.MolWt(mol)),
        "logp": float(crippen.MolLogP(mol)),
        "tpsa": float(rd_mol_descriptors.CalcTPSA(mol)),
        "hbd": float(lipinski.NumHDonors(mol)),
        "hba": float(lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": float(lipinski.NumRotatableBonds(mol)),
    }


def prepare_library(cfg: dict[str, Any], out_dir: Path, dry_run: bool = False, seed: int = 42) -> list[dict[str, Any]]:
    input_csv = Path(cfg.get("input_csv", ""))
    id_col = cfg.get("id_column", "compound_id")
    smiles_col = cfg.get("smiles_column", "smiles")
    label_col = cfg.get("label_column", "label")
    max_compounds = int(cfg.get("max_compounds", 50000))

    if not input_csv.exists():
        raise FileNotFoundError(f"Ligand library CSV not found: {input_csv}")

    compounds: list[dict[str, Any]] = []
    if dry_run:
        with input_csv.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if len(compounds) >= max_compounds:
                    break
                compounds.append(
                    {
                        "compound_id": row[id_col],
                        "smiles": row[smiles_col],
                        "label": int(row.get(label_col, 0)),
                        "mw": 300.0,
                        "logp": 2.5,
                        "tpsa": 80.0,
                        "hbd": 1.0,
                        "hba": 4.0,
                        "rotatable_bonds": 4.0,
                        "fingerprint": [0, 1] * 32,
                    }
                )
        return compounds

    chem, all_chem, crippen, descriptors, lipinski, rd_mol_descriptors = _safe_import_rdkit()

    sdf_path = out_dir / "library_prepared.sdf"
    writer = chem.SDWriter(str(sdf_path))

    with input_csv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if len(compounds) >= max_compounds:
                break
            smiles = row[smiles_col]
            mol = chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol = chem.AddHs(mol)
            embedding_result = all_chem.EmbedMolecule(mol, randomSeed=seed)
            if embedding_result != 0:
                continue
            all_chem.UFFOptimizeMolecule(mol)

            descriptor_values = _compute_descriptors(descriptors, crippen, lipinski, rd_mol_descriptors, mol)
            fingerprint = _fingerprint_bits(all_chem, mol)

            compound = {
                "compound_id": row[id_col],
                "smiles": smiles,
                "label": int(row.get(label_col, 0)),
                **descriptor_values,
                "fingerprint": fingerprint,
            }

            mol.SetProp("compound_id", compound["compound_id"])
            writer.write(mol)
            compounds.append(compound)

    writer.close()
    (out_dir / "library_prepared.json").write_text(json.dumps(compounds, indent=2), encoding="utf-8")
    return compounds
