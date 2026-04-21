from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


REQUIRED_TOP_LEVEL_KEYS = {
    "target",
    "pocket_detection",
    "library",
    "docking",
    "rescoring",
    "admet",
    "ranking",
}


def load_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path)
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Expected mapping in config file: {path}")
    missing = REQUIRED_TOP_LEVEL_KEYS - set(config.keys())
    if missing:
        missing_keys = ", ".join(sorted(missing))
        raise ValueError(f"Missing required config sections: {missing_keys}")
    return config
