# AlphaFold VS Pipeline (Publication-Grade SBVS)

A reproducible, modular, and explainable structure-based virtual screening (SBVS) platform for research use, MSc thesis workflows, and benchmark-driven method development.

## Features

- Real docking integration (AutoDock Vina/Smina CLI-compatible wrapper)
- RDKit ligand preparation (SMILES -> 3D conformers, descriptors, Morgan fingerprints)
- Pocket metadata pipeline (predefined boxes + fpocket integration path)
- Leakage-aware ML rescoring with LightGBM
- Benchmark metrics: ROC-AUC, PR-AUC, EF1%, calibration/Brier
- SHAP explainability artifacts (global + local outputs)
- ADMET profiling (Lipinski/Veber-style triage + toxicity alerts)
- Plotly + Streamlit research dashboard
- Deterministic runs (seeded execution + saved run config + JSON artifacts)

## Repository layout

- `src/alphafold_vs_pipeline/` - core pipeline and stage modules
- `configs/pipeline.yaml` - reproducible configuration
- `data/benchmark/demo_benchmark.csv` - demo benchmark dataset
- `dashboard/app.py` - publication-style Streamlit dashboard

## Installation

Base install (pipeline + tests):

```bash
pip install -e '.[dev]'
```

Dashboard support:

```bash
pip install -e '.[dev,dashboard]'
```

Full SBVS stack (RDKit/LightGBM/SHAP/Vina bindings):

```bash
pip install -e '.[dev,dashboard,full]'
```

## Run pipeline

```bash
alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --dry-run
```

Stage-wise execution:

```bash
alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --stage rescoring --dry-run
```

Real docking mode requires external binaries in PATH:

- `vina` or `smina`
- `obabel`
- optional `fpocket`

## Dashboard

```bash
streamlit run dashboard/app.py
```

Sections:

- Summary
- Docking Results
- ML Predictions
- SHAP Explainability
- ADMET Profile
- Ranking Table

## Reproducibility

Each run stores:

- `summary.json`
- `run_config.yaml` (serialized run configuration)
- `ml_metrics.json`
- `admet_profile.json`
- docking and model artifacts in output subdirectories

## Benchmarking

The default demo uses `data/benchmark/demo_benchmark.csv` with active/inactive labels for rescoring validation. Replace this with DUD-E or BindingDB subsets for publication benchmarking.

## Test

```bash
pytest -q
```
