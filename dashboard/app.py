from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
import streamlit as st


st.set_page_config(page_title="alphafold-vs-pipeline", layout="wide")
st.title("AlphaFold-guided virtual screening dashboard")


def discover_summary_paths(root: Path) -> list[str]:
    candidates = sorted(root.glob("outputs/**/summary.json"))
    return [str(p) for p in candidates if p.is_file()]


repo_root = Path.cwd()
detected_paths = discover_summary_paths(repo_root)

st.sidebar.subheader("Data source")
uploaded_file = st.sidebar.file_uploader("Upload summary JSON", type=["json"])
summary_path = st.sidebar.selectbox(
    "Detected run outputs",
    options=[""] + detected_paths,
    index=1 if detected_paths else 0,
    help="Choose a generated run output at outputs/<run>/summary.json",
)
manual_path = st.sidebar.text_input("Or enter summary JSON path", summary_path)

data: dict[str, Any]
if uploaded_file is not None:
    data = json.load(uploaded_file)
else:
    path = Path(manual_path.strip())
    if not manual_path.strip() or not path.exists():
        st.info(
            "Run the pipeline first and point this app to a real summary JSON, e.g. "
            "`alphafold-vs run --config configs/pipeline.yaml --output outputs/run_01` "
            "then use `outputs/run_01/summary.json`."
        )
        st.stop()
    data = json.loads(path.read_text(encoding="utf-8"))

hits = data.get("hits", [])
pockets = data.get("pockets", [])

st.subheader("Run summary")
col1, col2, col3, col4 = st.columns(4)
col1.metric("Target", str(data.get("target", "N/A")))
col2.metric("PDB ID", str(data.get("target_pdb_id", "N/A")))
col3.metric("Pocket count", len(pockets))
col4.metric("Hit count", len(hits))
st.caption(f"Run mode: {data.get('mode', 'N/A')}")

if pockets:
    st.subheader("Pocket summary")
    st.dataframe(pd.DataFrame(pockets), width="stretch", hide_index=True)

if hits:
    hits_df = pd.json_normalize(hits)
    st.subheader("Interactive hit database")
    sortable_columns = list(hits_df.columns)
    score_col = "docking_score" if "docking_score" in hits_df.columns else None
    rank_col = "ml_rescore" if "ml_rescore" in hits_df.columns else None

    filters = st.columns(3)
    top_n = filters[0].slider("Rows to display", min_value=1, max_value=len(hits_df), value=min(25, len(hits_df)))
    sort_col_default = score_col if score_col else sortable_columns[0]
    sort_col = filters[1].selectbox("Sort by", sortable_columns, index=sortable_columns.index(sort_col_default))
    ascending = filters[2].toggle("Ascending order", value=True)

    filtered_df = hits_df.sort_values(by=sort_col, ascending=ascending).head(top_n)
    st.dataframe(filtered_df, width="stretch", hide_index=True)

    fig = px.scatter_3d(
        filtered_df,
        x="docking_score",
        y="ml_rescore",
        z=[i + 1 for i in range(len(filtered_df))],
        hover_name="compound_id",
        title="Pose ranking (3D view)",
    )
    st.plotly_chart(fig, width="stretch")
    st.download_button("Download filtered hits (JSON)", data=filtered_df.to_json(orient="records", indent=2), file_name="hits.json")
