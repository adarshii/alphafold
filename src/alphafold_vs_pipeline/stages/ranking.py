from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def rank_and_export(compounds: list[dict[str, Any]], cfg: dict[str, Any], output_dir: Path) -> list[dict[str, Any]]:
    ranked = sorted(compounds, key=lambda x: x["ml_rescore"])
    ranked = ranked[: min(cfg.get("top_k", 100), len(ranked))]

    csv_path = output_dir / "hits.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["compound_id", "smiles", "docking_score", "ml_rescore"])
        writer.writeheader()
        for row in ranked:
            writer.writerow({key: row[key] for key in writer.fieldnames})

    if cfg.get("export_sdf", True):
        sdf_path = output_dir / "hits.sdf"
        with sdf_path.open("w", encoding="utf-8") as handle:
            for row in ranked:
                handle.write(f"{row['compound_id']}\n")
                handle.write("  alphafold-vs-pipeline\n\n")
                handle.write(">  <SMILES>\n")
                handle.write(f"{row['smiles']}\n\n")
                handle.write(">  <ML_RESCORE>\n")
                handle.write(f"{row['ml_rescore']}\n\n$$$$\n")
    return ranked
