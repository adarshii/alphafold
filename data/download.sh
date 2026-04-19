#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${ROOT_DIR}/raw"
mkdir -p "${RAW_DIR}"

echo "[1/3] Downloading PDB target 6LU7..."
curl -L "https://files.rcsb.org/download/6LU7.pdb" -o "${RAW_DIR}/6LU7.pdb"

echo "[2/3] Downloading ChEMBL actives (example query payload)..."
curl -L "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=CHEMBL3927&standard_type=IC50&standard_units=uM&limit=1000" \
  -o "${RAW_DIR}/chembl_mpro_ic50_uM.json"

echo "[3/3] Downloading ZINC15 drug-like decoy subset (placeholder URL)..."
wget -O "${RAW_DIR}/zinc15_druglike_10k.smi" \
  "https://files.docking.org/substances/subsets/drug-like/00/00.smi"

echo "Done. Files written under: ${RAW_DIR}"
