# alphafold-vs-pipeline

[![CI](https://github.com/adarshii/alphafold/actions/workflows/ci.yml/badge.svg)](https://github.com/adarshii/alphafold/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/adarshii/alphafold/branch/main/graph/badge.svg)](https://codecov.io/gh/adarshii/alphafold)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.TBD.svg)](https://doi.org/10.5281/zenodo.TBD)

End-to-end **AlphaFold-guided virtual screening** scaffold for drug discovery.

## Project goal

This repository provides a publication-ready baseline pipeline that:

1. Predicts or accepts a protein target structure (AlphaFold2 / ESMFold / ColabFold)
2. Detects druggable binding pockets (fpocket, DoGSiteScorer)
3. Prepares and filters compound libraries (RDKit, Open Babel, OpenEye OEChem, PyRx)
4. Runs batch molecular docking (AutoDock Vina, smina, Gnina)
5. Rescores poses with a LightGBM-ready feature interface (RF-Score + custom descriptors)
6. Filters compounds by ADMET-like criteria (Lipinski/PAINS hooks + external API hooks)
7. Ranks and exports hits with a Streamlit dashboard and CSV/SDF output

## Demo target and dataset

- **Target**: SARS-CoV-2 main protease (Mpro), PDB ID `6LU7`
- **Ligands**: ChEMBL Mpro actives (`IC50 < 10 µM`) + ZINC15 drug-like decoys (`10k`)
- Download commands are provided in:
  - `data/download.sh`

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .[dev]

# Run scaffolded pipeline
alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --dry-run

# Launch dashboard
streamlit run dashboard/app.py
```

## Repository layout

```text
configs/                     YAML pipeline configuration
dashboard/app.py             Streamlit + Plotly dashboard
data/download.sh             Demo data acquisition script
src/alphafold_vs_pipeline/   Core pipeline package
tests/                       Focused unit tests
```

## Citation

See `CITATION.cff`.
