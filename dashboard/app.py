from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


st.set_page_config(page_title="SBVS Research Dashboard", layout="wide")
st.title("Structure-Based Virtual Screening Dashboard")


def discover_summary_paths(root: Path) -> list[str]:
    candidates = sorted(root.glob("outputs/**/summary.json"))
    return [str(p) for p in candidates if p.is_file()]


def load_summary(uploaded: Any, path_text: str) -> dict[str, Any]:
    if uploaded is not None:
        raw = json.load(uploaded)
    else:
        path = Path(path_text.strip())
        if not path_text.strip() or not path.exists():
            raise FileNotFoundError("Choose a valid summary.json path or upload a summary file.")
        raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("summary.json must be a JSON object")
    raw.setdefault("hits", [])
    raw.setdefault("pockets", [])
    return raw


def _show_radar_for_hit(hit: dict[str, Any]) -> None:
    categories = ["mw", "logp", "tpsa", "hbd", "hba", "drug_likeness_score"]
    values = [float(hit.get(key, 0.0)) for key in categories]
    values.append(values[0])
    theta = categories + [categories[0]]
    fig = go.Figure(data=go.Scatterpolar(r=values, theta=theta, fill="toself"))
    fig.update_layout(template="plotly_white", title=f"ADMET Radar: {hit.get('compound_id', 'Top Hit')}")
    st.plotly_chart(fig, width="stretch")


repo_root = Path.cwd()
detected_paths = discover_summary_paths(repo_root)

st.sidebar.subheader("Data Source")
uploaded_file = st.sidebar.file_uploader("Upload summary JSON", type=["json"])
summary_path = st.sidebar.selectbox("Detected outputs", [""] + detected_paths, index=1 if detected_paths else 0)
manual_path = st.sidebar.text_input("Or enter summary path", summary_path)

try:
    data = load_summary(uploaded_file, manual_path)
except (FileNotFoundError, ValueError) as exc:
    st.info(str(exc))
    st.stop()

hits = data.get("hits", [])
hits_df = pd.json_normalize(hits) if hits else pd.DataFrame()
ml = data.get("ml", {})
metrics = ml.get("metrics", {})

st.header("Summary")
col1, col2, col3, col4 = st.columns(4)
col1.metric("Target", str(data.get("target", "N/A")))
col2.metric("PDB ID", str(data.get("target_pdb_id", "N/A")))
col3.metric("Pocket Count", len(data.get("pockets", [])))
col4.metric("Hit Count", len(hits))
st.caption(f"Mode: {data.get('mode', 'N/A')} | Seed: {data.get('seed', 'N/A')}")

st.header("Docking Results")
if hits_df.empty or "docking_score" not in hits_df.columns:
    st.info("No docking scores found.")
else:
    fig = px.histogram(hits_df, x="docking_score", nbins=20, title="Docking Score Distribution", template="plotly_white")
    st.plotly_chart(fig, width="stretch")

st.header("ML Predictions")
if hits_df.empty or "ml_rescore" not in hits_df.columns:
    st.info("No ML rescoring outputs found.")
else:
    fig = px.histogram(hits_df, x="ml_rescore", nbins=20, title="ML Rescore Distribution", template="plotly_white")
    st.plotly_chart(fig, width="stretch")

metric_cols = st.columns(4)
metric_cols[0].metric("ROC-AUC", f"{metrics.get('roc_auc', 0):.3f}" if metrics.get("roc_auc") is not None else "N/A")
metric_cols[1].metric("PR-AUC", f"{metrics.get('pr_auc', 0):.3f}" if metrics.get("pr_auc") is not None else "N/A")
metric_cols[2].metric("EF1%", f"{metrics.get('ef1', 0):.3f}" if metrics.get("ef1") is not None else "N/A")
metric_cols[3].metric("Brier", f"{metrics.get('brier', 0):.3f}" if metrics.get("brier") is not None else "N/A")

if isinstance(ml.get("calibration"), dict):
    cal = ml["calibration"]
    if cal.get("prob_pred") and cal.get("prob_true"):
        cal_df = pd.DataFrame({"predicted": cal["prob_pred"], "observed": cal["prob_true"]})
        fig = px.line(cal_df, x="predicted", y="observed", markers=True, title="Calibration Curve", template="plotly_white")
        st.plotly_chart(fig, width="stretch")

st.header("SHAP Explainability")
shap_meta = ml.get("shap", {}) if isinstance(ml, dict) else {}
summary_plot = Path(str(shap_meta.get("summary_plot", "")))
if summary_plot.exists():
    st.image(str(summary_plot), caption="SHAP Global Feature Importance")
else:
    st.info("SHAP summary plot not found for this run.")

st.header("ADMET Profile")
if hits_df.empty:
    st.info("No ADMET profile available.")
else:
    top_hit = hits_df.sort_values("ml_rescore", ascending=False).iloc[0].to_dict() if "ml_rescore" in hits_df.columns else hits_df.iloc[0].to_dict()
    _show_radar_for_hit(top_hit)

st.header("Ranking Table")
if hits_df.empty:
    st.info("No ranked hits available.")
else:
    sort_col = "ml_rescore" if "ml_rescore" in hits_df.columns else hits_df.columns[0]
    st.dataframe(hits_df.sort_values(sort_col, ascending=False), width="stretch", hide_index=True)
