from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def rank_and_export(compounds: list[dict[str, Any]], cfg: dict[str, Any], output_dir: Path) -> list[dict[str, Any]]:
    ranked = sorted(
        compounds,
        key=lambda x: (-float(x.get("ml_rescore", 0.0)), float(x.get("docking_score", 0.0))),
    )
    ranked = ranked[: min(int(cfg.get("top_k", 100)), len(ranked))]

    if cfg.get("export_csv", True):
        csv_path = output_dir / "hits.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            fieldnames = [
                "compound_id",
                "smiles",
                "docking_score",
                "ml_rescore",
                "drug_likeness_score",
                "admet_violations",
            ]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in ranked:
                writer.writerow({key: row.get(key, "") for key in fieldnames})

    if cfg.get("export_sdf", True):
        sdf_path = output_dir / "hits.sdf"
        with sdf_path.open("w", encoding="utf-8") as handle:
            for row in ranked:
                handle.write(f"{row['compound_id']}\n")
                handle.write("  alphafold-vs-pipeline\n\n")
                handle.write(">  <SMILES>\n")
                handle.write(f"{row['smiles']}\n\n")
                handle.write(">  <DOCKING_SCORE>\n")
                handle.write(f"{row.get('docking_score', '')}\n\n")
                handle.write(">  <ML_RESCORE>\n")
                handle.write(f"{row.get('ml_rescore', '')}\n\n")
                handle.write(">  <DRUG_LIKENESS_SCORE>\n")
                handle.write(f"{row.get('drug_likeness_score', '')}\n\n$$$$\n")

    return ranked
