import json
from pathlib import Path

from alphafold_vs_pipeline.cli import main


def test_cli_run_dry_run_outputs_files(tmp_path: Path) -> None:
    out_dir = tmp_path / "run"
    exit_code = main(["run", "--config", "configs/pipeline.yaml", "--output", str(out_dir), "--dry-run"])
    assert exit_code == 0
    summary = json.loads((out_dir / "summary.json").read_text(encoding="utf-8"))
    assert summary["target_pdb_id"] == "6LU7"
    assert (out_dir / "hits.csv").exists()
