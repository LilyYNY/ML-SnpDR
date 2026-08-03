"""Nested out-of-fold scoring for the canonical ML-SnpDR feature matrix."""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer


PROBABILITY_COLUMNS = [f"prob_C{i}" for i in range(1, 5)]
IDENTITY_COLUMNS = ["module_uid", "network", "method", "subtype", "module"]
QC_COLUMNS = ["feature_missing_count", "feature_qc_pass", "feature_qc_reason"]


@dataclass(frozen=True)
class ScoringConfig:
    random_state: int = 20260513
    outer_splits: int = 5
    outer_repeats: int = 5
    inner_splits: int = 3
    top_k: int = 10
    require_predicted_subtype_match: bool = True
    mode: str = "fast"
    n_jobs: int = 1

    def validate(self, class_counts: Mapping[str, int]) -> None:
        if self.mode not in {"fast", "paper"}:
            raise ValueError("mode must be 'fast' or 'paper'")
        for name in ("outer_splits", "outer_repeats", "inner_splits", "top_k"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive")
        smallest = min(class_counts.values())
        if smallest < self.outer_splits:
            raise ValueError("The smallest subtype has fewer rows than outer_splits")
        smallest_outer_train = smallest - int(np.ceil(smallest / self.outer_splits))
        if smallest_outer_train < self.inner_splits:
            raise ValueError("Outer training folds do not support inner_splits")


def _truthy(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes"})


def read_feature_matrix(feature_file: Path, schema_file: Path) -> tuple[pd.DataFrame, np.ndarray, List[str]]:
    features = pd.read_csv(feature_file, sep="\t")
    schema = json.loads(schema_file.read_text(encoding="utf-8"))
    feature_columns = list(schema["feature_columns"])
    required = IDENTITY_COLUMNS + QC_COLUMNS + feature_columns
    missing = [column for column in required if column not in features.columns]
    if missing:
        raise ValueError(f"module_features.tsv is missing: {', '.join(missing)}")
    if features["module_uid"].duplicated().any():
        raise ValueError("module_features.tsv contains duplicate module_uid values")
    if not _truthy(features["feature_qc_pass"]).all():
        raise ValueError("Every module must pass feature QC")
    matrix = features[feature_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not np.isfinite(matrix).all():
        raise ValueError("Feature matrix contains non-finite values")
    observed = sorted(features["subtype"].astype(str).unique())
    if observed != ["C1", "C2", "C3", "C4"]:
        raise ValueError(f"Expected subtypes C1-C4, found {observed}")
    return features, matrix, feature_columns


def parameter_grid(mode: str) -> List[Dict[str, Any]]:
    if mode == "fast":
        return [
            {
                "model__n_estimators": 220,
                "model__learning_rate": 0.08,
                "model__max_depth": 3,
                "model__subsample": 0.75,
                "model__min_samples_leaf": 1,
            }
        ]
    return [
        {
            "model__n_estimators": estimators,
            "model__learning_rate": rate,
            "model__max_depth": depth,
            "model__subsample": subsample,
            "model__min_samples_leaf": leaf,
        }
        for estimators in (120, 220)
        for rate in (0.025, 0.05, 0.08)
        for depth in (1, 2, 3)
        for subsample in (0.75, 1.0)
        for leaf in (1, 3)
    ]


def _pipeline(seed: int) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("model", GradientBoostingClassifier(random_state=seed)),
        ]
    )


def _plain_params(parameters: Mapping[str, Any]) -> Dict[str, Any]:
    return {key.replace("model__", "", 1): value for key, value in parameters.items()}


def _parameter_key(parameters: Mapping[str, Any]) -> str:
    return json.dumps(_plain_params(parameters), sort_keys=True, separators=(",", ":"))


def nested_oof_probabilities(
    matrix: np.ndarray,
    labels: np.ndarray,
    classes: Sequence[str],
    config: ScoringConfig,
) -> tuple[np.ndarray, List[Dict[str, Any]]]:
    splitter = RepeatedStratifiedKFold(
        n_splits=config.outer_splits,
        n_repeats=config.outer_repeats,
        random_state=config.random_state,
    )
    probability_sum = np.zeros((len(labels), len(classes)), dtype=float)
    prediction_count = np.zeros(len(labels), dtype=int)
    folds: List[Dict[str, Any]] = []
    grid = [
        {name: [value] for name, value in candidate.items()}
        for candidate in parameter_grid(config.mode)
    ]

    for fold_number, (train_index, test_index) in enumerate(splitter.split(matrix, labels), start=1):
        inner = StratifiedKFold(
            n_splits=config.inner_splits,
            shuffle=True,
            random_state=config.random_state + fold_number,
        )
        search = GridSearchCV(
            estimator=_pipeline(config.random_state + fold_number),
            param_grid=grid,
            scoring="roc_auc_ovr",
            cv=inner,
            n_jobs=config.n_jobs,
            refit=True,
            error_score="raise",
        )
        search.fit(matrix[train_index], labels[train_index])
        fold_probability = search.predict_proba(matrix[test_index])
        aligned = np.zeros((len(test_index), len(classes)), dtype=float)
        for local_index, class_name in enumerate(search.best_estimator_.classes_):
            aligned[:, list(classes).index(str(class_name))] = fold_probability[:, local_index]
        probability_sum[test_index] += aligned
        prediction_count[test_index] += 1
        folds.append(
            {
                "outer_fold": fold_number,
                "train_number": len(train_index),
                "test_number": len(test_index),
                "inner_best_score": float(search.best_score_),
                "best_params": _plain_params(search.best_params_),
            }
        )

    if np.any(prediction_count != config.outer_repeats):
        raise RuntimeError("Repeated OOF coverage is incomplete")
    return probability_sum / prediction_count[:, None], folds


def build_score_tables(
    features: pd.DataFrame,
    probabilities: np.ndarray,
    classes: Sequence[str],
    config: ScoringConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    class_index = {name: index for index, name in enumerate(classes)}
    true_subtype = features["subtype"].astype(str).to_numpy()
    predicted = np.asarray(classes)[np.argmax(probabilities, axis=1)]
    target_probability = np.array(
        [probabilities[row, class_index[subtype]] for row, subtype in enumerate(true_subtype)]
    )
    competing_probability = np.array(
        [
            np.max(np.delete(probabilities[row], class_index[subtype]))
            for row, subtype in enumerate(true_subtype)
        ]
    )
    scores = features[IDENTITY_COLUMNS].copy()
    scores["true_subtype"] = true_subtype
    scores["predicted_subtype"] = predicted
    for index, column in enumerate(PROBABILITY_COLUMNS):
        scores[column] = probabilities[:, index]
    scores["target_subtype_probability"] = target_probability
    scores["probability_margin"] = target_probability - competing_probability
    scores["rank_in_subtype"] = (
        scores.sort_values(
            ["true_subtype", "target_subtype_probability", "probability_margin", "module_uid"],
            ascending=[True, False, False, True],
        )
        .groupby("true_subtype")
        .cumcount()
        .add(1)
        .sort_index()
        .astype(int)
    )
    scores["score_type"] = "nested_oof_probability"
    scores["model_version"] = "gradient_boosting_nested_oof_v1"

    eligible = scores.copy()
    if config.require_predicted_subtype_match:
        eligible = eligible[eligible["predicted_subtype"] == eligible["true_subtype"]]
    top_parts = []
    for subtype in classes:
        group = eligible[eligible["true_subtype"] == subtype].sort_values(
            ["rank_in_subtype", "module_uid"]
        )
        if len(group) < config.top_k:
            raise ValueError(
                f"Subtype {subtype} has {len(group)} eligible modules; top_k={config.top_k}"
            )
        group = group.head(config.top_k).copy()
        group["target_subtype"] = subtype
        group["top_k"] = config.top_k
        group["ml_gate"] = True
        group["ml_gate_reason"] = "passed_configured_ml_gates"
        top_parts.append(group)
    return scores, pd.concat(top_parts, ignore_index=True)


def score_metrics(labels: np.ndarray, probabilities: np.ndarray, classes: Sequence[str]) -> Dict[str, float]:
    predicted = np.asarray(classes)[np.argmax(probabilities, axis=1)]
    return {
        "macro_ovr_auc": float(
            roc_auc_score(labels, probabilities, labels=list(classes), multi_class="ovr", average="macro")
        ),
        "accuracy": float(accuracy_score(labels, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predicted)),
        "macro_f1": float(f1_score(labels, predicted, average="macro")),
    }


def run_nested_scoring(
    feature_file: Path,
    schema_file: Path,
    output_dir: Path,
    config: ScoringConfig,
) -> Dict[str, Path]:
    feature_file = feature_file.resolve(strict=True)
    schema_file = schema_file.resolve(strict=True)
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    features, matrix, feature_columns = read_feature_matrix(feature_file, schema_file)
    labels = features["subtype"].astype(str).to_numpy()
    classes = sorted(np.unique(labels).tolist())
    class_counts = Counter(labels)
    config.validate(class_counts)
    probabilities, folds = nested_oof_probabilities(matrix, labels, classes, config)
    scores, top = build_score_tables(features, probabilities, classes, config)
    metrics = score_metrics(labels, probabilities, classes)

    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}-", dir=output_dir.parent))
    try:
        scores.to_csv(staging / "ml_scores.tsv", sep="\t", index=False)
        top.to_csv(staging / f"ml_top{config.top_k}.tsv", sep="\t", index=False)
        pd.DataFrame(
            [
                {
                    **{key: value for key, value in fold.items() if key != "best_params"},
                    "best_params": json.dumps(fold["best_params"], sort_keys=True),
                }
                for fold in folds
            ]
        ).to_csv(staging / "nested_oof_folds.tsv", sep="\t", index=False)

        chosen_key = Counter(_parameter_key(fold["best_params"]) for fold in folds).most_common(1)[0][0]
        chosen = json.loads(chosen_key)
        final_model = _pipeline(config.random_state)
        final_model.set_params(**{f"model__{key}": value for key, value in chosen.items()})
        final_model.fit(matrix, labels)
        joblib.dump(final_model, staging / "fitted_model.joblib")
        metadata = {
            "config": asdict(config),
            "classes": classes,
            "class_counts": dict(sorted(class_counts.items())),
            "feature_columns": feature_columns,
            "n_modules": len(features),
            "n_features": len(feature_columns),
            "metrics": metrics,
            "final_model_params": chosen,
            "score_type": "nested_oof_probability",
        }
        (staging / "ml_scoring_metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
        )
        staging.replace(output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {
        "scores_file": output_dir / "ml_scores.tsv",
        "top_file": output_dir / f"ml_top{config.top_k}.tsv",
        "metadata_file": output_dir / "ml_scoring_metadata.json",
        "model_file": output_dir / "fitted_model.joblib",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, type=Path, help="Step-6A module_features.tsv")
    parser.add_argument("--schema", required=True, type=Path, help="Step-6A feature_schema.json")
    parser.add_argument("--output", required=True, type=Path, help="New step-6B output directory")
    parser.add_argument("--mode", choices=("fast", "paper"), default="fast")
    parser.add_argument("--seed", type=int, default=20260513)
    parser.add_argument("--outer-splits", type=int, default=5)
    parser.add_argument("--outer-repeats", type=int, default=5)
    parser.add_argument("--inner-splits", type=int, default=3)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument(
        "--allow-predicted-subtype-mismatch",
        action="store_true",
        help="Allow a module into Top-K when its predicted class differs from its true subtype",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    config = ScoringConfig(
        random_state=args.seed,
        outer_splits=args.outer_splits,
        outer_repeats=args.outer_repeats,
        inner_splits=args.inner_splits,
        top_k=args.top_k,
        require_predicted_subtype_match=not args.allow_predicted_subtype_mismatch,
        mode=args.mode,
        n_jobs=args.n_jobs,
    )
    outputs = run_nested_scoring(args.features, args.schema, args.output, config)
    print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
