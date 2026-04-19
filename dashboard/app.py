from __future__ import annotations

import json
from pathlib import Path

import plotly.express as px
import streamlit as st


st.set_page_config(page_title="alphafold-vs-pipeline", layout="wide")
st.title("AlphaFold-guided virtual screening dashboard")

summary_path = st.sidebar.text_input("Summary JSON path", "outputs/demo/summary.json")
path = Path(summary_path)

if not path.exists():
    st.info("Run `alphafold-vs run --config configs/pipeline.yaml --output outputs/demo --dry-run` first.")
    st.stop()

data = json.loads(path.read_text(encoding="utf-8"))
hits = data.get("hits", [])

st.subheader("Run summary")
st.json(
    {
        "target": data.get("target"),
        "target_pdb_id": data.get("target_pdb_id"),
        "pocket_count": len(data.get("pockets", [])),
        "hit_count": len(hits),
        "mode": data.get("mode"),
    }
)

if hits:
    fig = px.scatter_3d(
        hits,
        x="docking_score",
        y="ml_rescore",
        z=[i + 1 for i in range(len(hits))],
        hover_name="compound_id",
        title="Pose ranking (mock 3D view)",
    )
    st.plotly_chart(fig, use_container_width=True)
    st.download_button("Download hits (JSON)", data=json.dumps(hits, indent=2), file_name="hits.json")
