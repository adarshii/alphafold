from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
from sklearn.calibration import calibration_curve
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split


def _safe_import_lightgbm() -> Any:
    try:
        from lightgbm import LGBMClassifier
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("LightGBM is required for rescoring; install with full extras.") from exc
    return LGBMClassifier


def _safe_import_shap() -> Any:
    try:
        import shap
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("SHAP is required when perform_shap is enabled.") from exc
    return shap


def _feature_columns(compounds: list[dict[str, Any]]) -> list[str]:
    candidates = ["docking_score", "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds"]
    return [col for col in candidates if col in compounds[0]]


def _to_matrix(compounds: list[dict[str, Any]], columns: list[str]) -> np.ndarray:
    return np.asarray([[float(c[col]) for col in columns] for c in compounds], dtype=float)


def _enrichment_factor_at_1_percent(y_true: np.ndarray, y_score: np.ndarray) -> float:
    n_total = len(y_true)
    n_top = max(1, int(np.ceil(0.01 * n_total)))
    order = np.argsort(-y_score)
    top_true = y_true[order[:n_top]]
    observed = float(np.sum(top_true)) / n_top
    baseline = float(np.sum(y_true)) / n_total if n_total else 0.0
    if baseline == 0.0:
        return 0.0
    return observed / baseline


def rescore_poses(
    compounds: list[dict[str, Any]],
    cfg: dict[str, Any],
    benchmark_cfg: dict[str, Any],
    out_dir: Path,
    dry_run: bool = False,
    seed: int = 42,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not compounds:
        return compounds, {}

    model_dir = out_dir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)

    columns = _feature_columns(compounds)
    X = _to_matrix(compounds, columns)
    y = np.asarray([int(c.get("label", 0)) for c in compounds], dtype=int)

    if dry_run or len(np.unique(y)) < 2:
        probabilities = 1 / (1 + np.exp(X[:, 0]))
        for item, score in zip(compounds, probabilities, strict=True):
            item["ml_rescore"] = float(round(score, 6))
        artifacts = {
            "features": columns,
            "metrics": {"roc_auc": None, "pr_auc": None, "ef1": None, "brier": None},
            "benchmark_dataset": benchmark_cfg.get("dataset_name", "unlabeled"),
        }
        (out_dir / "ml_metrics.json").write_text(json.dumps(artifacts, indent=2), encoding="utf-8")
        return compounds, artifacts

    lgbm_classifier = _safe_import_lightgbm()
    class_weight = cfg.get("class_weight", "balanced")
    model = lgbm_classifier(random_state=seed, class_weight=class_weight, n_estimators=200, learning_rate=0.05)

    test_size = max(0.1, 1 - float(cfg.get("train_fraction", 0.8)))
    X_train, X_test, y_train, y_test, train_idx, test_idx = train_test_split(
        X,
        y,
        np.arange(len(compounds)),
        test_size=test_size,
        random_state=seed,
        stratify=y,
    )

    cv_folds = int(cfg.get("cv_folds", 5))
    cv = StratifiedKFold(n_splits=min(cv_folds, np.bincount(y_train).min()), shuffle=True, random_state=seed)
    cv_scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc")

    model.fit(X_train, y_train)
    full_probs = model.predict_proba(X)[:, 1]
    test_probs = model.predict_proba(X_test)[:, 1]

    for idx, score in enumerate(full_probs):
        compounds[idx]["ml_rescore"] = float(score)

    roc_auc = float(roc_auc_score(y_test, test_probs))
    pr_auc = float(average_precision_score(y_test, test_probs))
    ef1 = float(_enrichment_factor_at_1_percent(y_test, test_probs))
    brier = float(brier_score_loss(y_test, test_probs))

    prob_true, prob_pred = calibration_curve(y_test, test_probs, n_bins=10)
    calibration = {"prob_true": prob_true.tolist(), "prob_pred": prob_pred.tolist()}

    model_path = Path(cfg.get("save_model_path", model_dir / "lightgbm_rescorer.joblib"))
    if not model_path.is_absolute():
        model_path = out_dir / model_path
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, model_path)

    artifacts: dict[str, Any] = {
        "features": columns,
        "cv_roc_auc_mean": float(np.mean(cv_scores)),
        "cv_roc_auc_std": float(np.std(cv_scores)),
        "metrics": {"roc_auc": roc_auc, "pr_auc": pr_auc, "ef1": ef1, "brier": brier},
        "calibration": calibration,
        "benchmark_dataset": benchmark_cfg.get("dataset_name", "custom"),
        "train_size": int(len(train_idx)),
        "test_size": int(len(test_idx)),
        "model_path": str(model_path),
    }

    if cfg.get("perform_shap", True):
        artifacts["shap"] = _generate_shap_outputs(model, X_train, columns, compounds, out_dir)

    (out_dir / "ml_metrics.json").write_text(json.dumps(artifacts, indent=2), encoding="utf-8")
    return compounds, artifacts


def _generate_shap_outputs(
    model: Any,
    X_train: np.ndarray,
    columns: list[str],
    compounds: list[dict[str, Any]],
    out_dir: Path,
) -> dict[str, Any]:
    shap = _safe_import_shap()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shap_dir = out_dir / "shap"
    shap_dir.mkdir(parents=True, exist_ok=True)

    explainer = shap.TreeExplainer(model)
    sample = X_train[: min(200, len(X_train))]
    shap_values = explainer.shap_values(sample)

    summary_png = shap_dir / "shap_summary.png"
    plt.figure(figsize=(8, 5))
    shap.summary_plot(shap_values, sample, feature_names=columns, show=False)
    plt.tight_layout()
    plt.savefig(summary_png, dpi=200)
    plt.close()

    top_hit = max(compounds, key=lambda item: float(item.get("ml_rescore", 0.0)))
    local_json = shap_dir / "local_explanation_top_hit.json"
    local_json.write_text(json.dumps({"compound_id": top_hit["compound_id"], "features": columns}, indent=2), encoding="utf-8")

    return {
        "summary_plot": str(summary_png),
        "local_explanation": str(local_json),
    }
