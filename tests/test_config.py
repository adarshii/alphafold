from pathlib import Path

from alphafold_vs_pipeline.config import load_config


def test_load_config_reads_target_block() -> None:
    config = load_config(Path("configs/pipeline.yaml"))
    assert config["target"]["pdb_id"] == "6LU7"
    assert "fpocket" in config["pocket_detection"]["tools"]
