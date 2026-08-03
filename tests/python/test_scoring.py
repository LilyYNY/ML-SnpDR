import json

import numpy as np
import pandas as pd

from mlsnpdr.scoring import ScoringConfig, run_nested_scoring


def test_nested_scoring_writes_all_module_scores_and_top_k(tmp_path):
    rows = []
    rng = np.random.default_rng(20260513)
    for class_index, subtype in enumerate(("C1", "C2", "C3", "C4"), start=1):
        for module_index in range(1, 9):
            rows.append(
                {
                    "module_uid": f"string__louvain__{subtype}__M{module_index}",
                    "network": "String",
                    "method": "Louvain",
                    "subtype": subtype,
                    "module": f"M{module_index}",
                    "feature_1": class_index + rng.normal(0, 0.05),
                    "feature_2": class_index**2 + rng.normal(0, 0.05),
                    "feature_3": (5 - class_index) + rng.normal(0, 0.05),
                    "feature_missing_count": 0,
                    "feature_qc_pass": True,
                    "feature_qc_reason": "pass",
                }
            )
    features = pd.DataFrame(rows)
    feature_file = tmp_path / "module_features.tsv"
    schema_file = tmp_path / "feature_schema.json"
    output_dir = tmp_path / "scores"
    features.to_csv(feature_file, sep="\t", index=False)
    schema_file.write_text(
        json.dumps({"feature_columns": ["feature_1", "feature_2", "feature_3"]}),
        encoding="utf-8",
    )

    outputs = run_nested_scoring(
        feature_file,
        schema_file,
        output_dir,
        ScoringConfig(
            outer_splits=2,
            outer_repeats=1,
            inner_splits=2,
            top_k=2,
            require_predicted_subtype_match=False,
            mode="fast",
            n_jobs=1,
        ),
    )

    scores = pd.read_csv(outputs["scores_file"], sep="\t")
    top = pd.read_csv(outputs["top_file"], sep="\t")
    metadata = json.loads(outputs["metadata_file"].read_text(encoding="utf-8"))
    assert len(scores) == 32
    assert len(top) == 8
    assert top.groupby("target_subtype").size().to_dict() == {
        "C1": 2,
        "C2": 2,
        "C3": 2,
        "C4": 2,
    }
    assert np.allclose(scores[[f"prob_C{i}" for i in range(1, 5)]].sum(axis=1), 1)
    assert metadata["n_features"] == 3
    assert outputs["model_file"].exists()
