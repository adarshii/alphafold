from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
import streamlit as st


st.set_page_config(page_title="alphafold-vs-pipeline", layout="wide")
st.title("AlphaFold-guided virtual screening dashboard")
PLACEHOLDER_VALUE = "N/A"


def sort_hits_dataframe(df: pd.DataFrame, sort_col: str, ascending: bool) -> pd.DataFrame:
    """Sort hit rows and fallback to string ordering for mixed/non-orderable values.

    If `sort_values` raises ``TypeError`` (for example, when a sort column mixes
    incomparable scalar types), values are converted into a temporary string key
    while preserving missing values (`None`, `NaN`, `NaT`) as empty strings.
    The helper key column is removed before returning the sorted dataframe.
    """
    try:
        return df.sort_values(by=sort_col, ascending=ascending)
    except TypeError:
        return df.assign(
            _sort_key=df[sort_col].map(lambda value: "" if pd.isna(value) else str(value))
        ).sort_values(by="_sort_key", ascending=ascending).drop(columns="_sort_key")


def discover_summary_paths(root: Path) -> list[str]:
    candidates = sorted(root.glob("outputs/**/summary.json"))
    return [str(p) for p in candidates if p.is_file()]


def normalize_summary_data(raw_data: Any) -> dict[str, Any]:
    """Normalize uploaded/loaded JSON into a dashboard-compatible summary dict.

    Accepts either:
    - a summary object (dict), where missing/non-list `hits` and `pockets` are coerced to lists
    - a list of hit objects, wrapped into a minimal summary payload

    Returns:
        dict[str, Any]: Normalized summary payload with at least `target`, `target_pdb_id`,
        `pockets`, and `hits`; includes `mode` for list-based uploaded hit payloads.

    Raises:
        ValueError: If the JSON shape is neither a summary object nor a list of hit objects.
    """
    if isinstance(raw_data, dict):
        normalized = dict(raw_data)
        if isinstance(normalized.get("hits"), list):
            normalized["hits"] = [item for item in normalized["hits"] if isinstance(item, dict)]
        else:
            normalized["hits"] = []
        if not isinstance(normalized.get("pockets"), list):
            normalized["pockets"] = []
        return normalized
    # Intentionally validate every row so downstream table/plot code only sees dict-like hit records.
    if isinstance(raw_data, list) and all(isinstance(item, dict) for item in raw_data):
        return {
            "target": PLACEHOLDER_VALUE,
            "target_pdb_id": PLACEHOLDER_VALUE,
            # Keep structure key present to match summary.json shape used across the app/docs.
            "structure": {},
            "pockets": [],
            "hits": raw_data,
            "mode": "uploaded-hits-only",
        }
    raise ValueError("Unsupported JSON format. Upload a summary object or a list of hit objects.")


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
raw_data: Any
if uploaded_file is not None:
    raw_data = json.load(uploaded_file)
else:
    path = Path(manual_path.strip())
    if not manual_path.strip() or not path.exists():
        st.info(
            "Run the pipeline first and point this app to a real summary JSON, e.g. "
            "`alphafold-vs run --config configs/pipeline.yaml --output outputs/run_01` "
            "then use `outputs/run_01/summary.json`."
        )
        st.stop()
    raw_data = json.loads(path.read_text(encoding="utf-8"))

try:
    data = normalize_summary_data(raw_data)
except ValueError as exc:
    st.error(str(exc))
    st.stop()

hits = data.get("hits", [])
pockets = data.get("pockets", [])

st.subheader("Run summary")
col1, col2, col3, col4 = st.columns(4)
col1.metric("Target", str(data.get("target", PLACEHOLDER_VALUE)))
col2.metric("PDB ID", str(data.get("target_pdb_id", PLACEHOLDER_VALUE)))
col3.metric("Pocket count", len(pockets))
col4.metric("Hit count", len(hits))
st.caption(f"Run mode: {data.get('mode', PLACEHOLDER_VALUE)}")

if pockets:
    st.subheader("Pocket summary")
    st.dataframe(pd.DataFrame(pockets), width="stretch", hide_index=True)

if hits:
    hits_df = pd.json_normalize(hits)
    st.subheader("Interactive hit database")
    if hits_df.empty:
        st.info("No tabular hit rows are available to display.")
    else:
        sortable_columns = list(hits_df.columns)
        score_col = "docking_score" if "docking_score" in hits_df.columns else None
        rank_col = "ml_rescore" if "ml_rescore" in hits_df.columns else None

        filters = st.columns(3)
        top_n = filters[0].slider("Rows to display", min_value=1, max_value=len(hits_df), value=min(25, len(hits_df)))
        sort_col_default = score_col if score_col else sortable_columns[0]
        sort_col = filters[1].selectbox("Sort by", sortable_columns, index=sortable_columns.index(sort_col_default))
        default_ascending = sort_col in {"docking_score", "ml_rescore"}
        ascending = filters[2].toggle("Sort ascending", value=default_ascending)

        filtered_df = sort_hits_dataframe(hits_df, sort_col, ascending).head(top_n)
        st.dataframe(filtered_df, width="stretch", hide_index=True)

        if {"docking_score", "ml_rescore"}.issubset(filtered_df.columns):
            fig = px.scatter_3d(
                filtered_df,
                x="docking_score",
                y="ml_rescore",
                z=list(range(1, len(filtered_df) + 1)),
                hover_name="compound_id",
                title="Pose ranking (3D view)",
            )
            fig.update_layout(autosize=True)
            st.plotly_chart(fig)
        else:
            st.info("3D ranking view requires `docking_score` and `ml_rescore` columns.")
        st.download_button(
            "Download filtered hits (JSON)",
            data=json.dumps(filtered_df.to_dict(orient="records"), indent=2),
            file_name="hits.json",
        )
