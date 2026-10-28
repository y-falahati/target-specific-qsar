
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)


# --------------------------------------------------
# Configuration
# --------------------------------------------------
RANDOM_SEED = 42
N_ITER = 20

INPUT_PATH = Path("data/processed/egfr_scaffold_split.csv")
BASELINE_METRICS_PATH = Path(
    "results/tables/egfr_baseline_validation_metrics.csv"
)

TABLES_DIR = Path("results/tables")
MODELS_DIR = Path("models")

FEATURES = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotatableBonds",
    "RingCount",
    "HeavyAtomCount",
    "FractionCSP3",
    "AromaticRings",
    "AliphaticRings",
    "HeteroatomCount",
    "FormalCharge",
]

TARGET = "pIC50"


def calculate_metrics(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "MAE": mean_absolute_error(y_true, y_pred),
    }


def main():
    df = pd.read_csv(INPUT_PATH)

    required = FEATURES + [
        TARGET,
        "split",
        "scaffold",
        "molecule_chembl_id",
    ]

    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    if df[required].isna().any().any():
        raise ValueError("Missing values found in modeling columns.")

    train = df[df["split"] == "train"].copy()
    validation = df[df["split"] == "validation"].copy()

    if train.empty or validation.empty:
        raise ValueError("Train or validation set is empty.")

    # Confirm that each scaffold belongs to only one outer split.
    split_counts = df.groupby("scaffold")["split"].nunique()
    if (split_counts > 1).any():
        raise ValueError("Scaffold leakage detected in outer split.")

    X_train = train[FEATURES]
    y_train = train[TARGET]
    groups_train = train["scaffold"]

    X_val = validation[FEATURES]
    y_val = validation[TARGET]

    # Grouped CV keeps the same scaffold out of both sides
    # of each internal train/validation fold.
    cv = GroupKFold(n_splits=5)

    parameter_space = {
        "n_estimators": [200, 300, 500, 700],
        "max_depth": [None, 10, 20, 30],
        "min_samples_split": [2, 5, 10],
        "min_samples_leaf": [1, 2, 4, 8],
        "max_features": [1.0, "sqrt", 0.7],
        "bootstrap": [True],
    }

    search = RandomizedSearchCV(
        estimator=RandomForestRegressor(
            random_state=RANDOM_SEED,
            n_jobs=1,
        ),
        param_distributions=parameter_space,
        n_iter=N_ITER,
        scoring="neg_root_mean_squared_error",
        cv=cv,
        random_state=RANDOM_SEED,
        n_jobs=-1,
        verbose=2,
        refit=True,
        return_train_score=True,
    )

    print("=" * 65)
    print("RANDOM FOREST HYPERPARAMETER TUNING")
    print("=" * 65)
    print(f"Training molecules: {len(train):,}")
    print(f"Validation molecules: {len(validation):,}")
    print(f"Unique training scaffolds: {groups_train.nunique():,}")
    print("Search: 20 randomized configurations, 5-fold GroupKFold")
    print("Test set will not be used.")

    # Search and refit using training data only.
    search.fit(X_train, y_train, groups=groups_train)

    best_model = search.best_estimator_

    # Evaluate the selected model once on the outer validation set.
    val_predictions = best_model.predict(X_val)
    tuned_metrics = calculate_metrics(y_val, val_predictions)

    print("\nBest parameters:")
    print(search.best_params_)

    print(
        "\nBest grouped-CV RMSE: "
        f"{-search.best_score_:.4f}"
    )

    print("\nTuned model — outer validation:")
    for metric, value in tuned_metrics.items():
        print(f"{metric}: {value:.4f}")

    # Compare against the previously saved baseline result.
    baseline = pd.read_csv(BASELINE_METRICS_PATH)
    baseline_rf = baseline.loc[
        baseline["model"] == "RandomForest"
    ]

    if baseline_rf.empty:
        raise ValueError("Random Forest baseline metrics not found.")

    baseline_row = baseline_rf.iloc[0]

    comparison = pd.DataFrame([
        {
            "model": "RandomForest_baseline",
            "R2": baseline_row["R2"],
            "RMSE": baseline_row["RMSE"],
            "MAE": baseline_row["MAE"],
        },
        {
            "model": "RandomForest_tuned",
            **tuned_metrics,
        },
    ])

    # Save tuning results and the fitted candidate model.
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    cv_results = pd.DataFrame(search.cv_results_)
    cv_results.to_csv(
        TABLES_DIR / "egfr_rf_tuning_cv_results.csv",
        index=False,
    )

    comparison.to_csv(
        TABLES_DIR / "egfr_rf_tuning_comparison.csv",
        index=False,
    )

    prediction_output = validation[
        ["molecule_chembl_id", TARGET]
    ].copy()

    prediction_output = prediction_output.rename(
        columns={TARGET: "y_true_pIC50"}
    )
    prediction_output["y_pred_pIC50"] = val_predictions
    prediction_output.to_csv(
        TABLES_DIR / "egfr_rf_tuned_validation_predictions.csv",
        index=False,
    )

    model_path = MODELS_DIR / "randomforest_tuned_candidate.joblib"
    joblib.dump(best_model, model_path)

    print("\nBaseline vs tuned validation performance:")
    print(comparison.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    ))

    print("\nSaved:")
    print(TABLES_DIR / "egfr_rf_tuning_cv_results.csv")
    print(TABLES_DIR / "egfr_rf_tuning_comparison.csv")
    print(TABLES_DIR / "egfr_rf_tuned_validation_predictions.csv")
    print(model_path)

    print("\nTest set was not used.")
    print("QC status: PASS")


if __name__ == "__main__":
    main()
