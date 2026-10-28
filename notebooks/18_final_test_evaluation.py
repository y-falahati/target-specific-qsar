
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)


# --------------------------------------------------
# 1. Project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

DATA_PATH = ROOT / "data" / "processed" / "egfr_scaffold_split.csv"
MODEL_DIR = ROOT / "models"
RESULTS_DIR = ROOT / "results" / "tables"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 2. Fixed evaluation settings
# --------------------------------------------------

TARGET = "pIC50"

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

MODELS = {
    "RandomForest_Baseline": (
        MODEL_DIR / "randomforest_baseline.joblib"
    ),
    "RandomForest_Tuned": (
        MODEL_DIR / "randomforest_tuned_candidate.joblib"
    ),
}


# --------------------------------------------------
# 3. Load the previously created split dataset
# --------------------------------------------------

df = pd.read_csv(DATA_PATH)

# Detect the split column without changing the data.
split_candidates = ["split", "dataset_split", "set", "subset"]
split_column = next(
    (col for col in split_candidates if col in df.columns),
    None,
)

if split_column is None:
    raise ValueError(
        "Could not identify the split column. "
        f"Available columns: {df.columns.tolist()}"
    )

test_mask = (
    df[split_column]
    .astype(str)
    .str.strip()
    .str.lower()
    .eq("test")
)

test_df = df.loc[test_mask].copy()

if test_df.empty:
    raise ValueError(
        "No Test rows were found. Check the split column and labels."
    )

required_columns = FEATURES + [TARGET]

missing_columns = [
    col for col in required_columns
    if col not in test_df.columns
]

if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}"
    )

if test_df[required_columns].isna().any().any():
    raise ValueError(
        "Missing values found in Test features or target. "
        "Investigate data quality before evaluation."
    )

X_test = test_df[FEATURES]
y_test = test_df[TARGET].astype(float)

if not np.isfinite(X_test.to_numpy(dtype=float)).all():
    raise ValueError("Non-finite values found in Test features.")

if not np.isfinite(y_test.to_numpy()).all():
    raise ValueError("Non-finite values found in Test target.")

print("=" * 60)
print("FINAL TEST EVALUATION")
print("=" * 60)
print(f"Dataset: {DATA_PATH}")
print(f"Split column: {split_column}")
print(f"Number of Test compounds: {len(test_df)}")
print(f"Number of features: {len(FEATURES)}")
print("No model fitting or hyperparameter tuning will be performed.")


# --------------------------------------------------
# 4. Evaluate saved models without retraining
# --------------------------------------------------

metrics_rows = []
prediction_frames = []

for model_name, model_path in MODELS.items():

    if not model_path.exists():
        raise FileNotFoundError(
            f"Saved model not found: {model_path}"
        )

    model = joblib.load(model_path)

    # Confirm compatibility with the saved feature set.
    if hasattr(model, "n_features_in_"):
        if model.n_features_in_ != len(FEATURES):
            raise ValueError(
                f"{model_name}: expected "
                f"{model.n_features_in_} features, "
                f"but FEATURES contains {len(FEATURES)}."
            )

    if hasattr(model, "feature_names_in_"):
        saved_features = list(model.feature_names_in_)
        if saved_features != FEATURES:
            raise ValueError(
                f"{model_name}: feature names/order differ "
                "from the configured FEATURES list.\n"
                f"Saved: {saved_features}\n"
                f"Configured: {FEATURES}"
            )

    # IMPORTANT: predict only. Do not call fit() or refit().
    y_pred = np.asarray(model.predict(X_test), dtype=float)

    if len(y_pred) != len(y_test):
        raise ValueError(
            f"{model_name}: prediction count does not match Test rows."
        )

    if not np.isfinite(y_pred).all():
        raise ValueError(
            f"{model_name}: non-finite predictions detected."
        )

    residuals = y_test.to_numpy() - y_pred
    absolute_errors = np.abs(residuals)

    metrics_rows.append({
        "model": model_name,
        "n_test": len(y_test),
        "R2": r2_score(y_test, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_test, y_pred)),
        "MAE": mean_absolute_error(y_test, y_pred),
        "mean_residual_observed_minus_predicted": residuals.mean(),
        "median_absolute_error": np.median(absolute_errors),
        "max_absolute_error": absolute_errors.max(),
    })

    prediction_frame = pd.DataFrame({
        "molecule_chembl_id": (
            test_df["molecule_chembl_id"].to_numpy()
            if "molecule_chembl_id" in test_df.columns
            else test_df.index.to_numpy()
        ),
        "observed_pIC50": y_test.to_numpy(),
        f"{model_name}_predicted_pIC50": y_pred,
        f"{model_name}_residual": residuals,
        f"{model_name}_absolute_error": absolute_errors,
    })

    prediction_frames.append(prediction_frame)

    print(f"\n{model_name}")
    print("-" * 40)
    print(f"R²:   {r2_score(y_test, y_pred):.4f}")
    print(
        f"RMSE: "
        f"{np.sqrt(mean_squared_error(y_test, y_pred)):.4f}"
    )
    print(f"MAE:  {mean_absolute_error(y_test, y_pred):.4f}")
    print(
        "Mean residual (observed - predicted): "
        f"{residuals.mean():.4f}"
    )


# --------------------------------------------------
# 5. Save final metrics and compound-level predictions
# --------------------------------------------------

metrics_df = pd.DataFrame(metrics_rows)

metrics_path = RESULTS_DIR / "egfr_final_test_metrics.csv"
metrics_df.to_csv(metrics_path, index=False)

predictions_df = prediction_frames[0]

for frame in prediction_frames[1:]:
    predictions_df = predictions_df.merge(
        frame,
        on=["molecule_chembl_id", "observed_pIC50"],
        how="inner",
        validate="one_to_one",
    )

predictions_path = (
    RESULTS_DIR / "egfr_final_test_predictions.csv"
)
predictions_df.to_csv(predictions_path, index=False)


# --------------------------------------------------
# 6. Report the final comparison
# --------------------------------------------------

print("\n" + "=" * 60)
print("FINAL TEST METRICS")
print("=" * 60)
print(metrics_df.to_string(index=False))

print("\nSaved files:")
print(metrics_path)
print(predictions_path)

print("\nEvaluation completed.")
print(
    "The Test set was used only for final prediction and reporting. "
    "Do not use these results to tune or select another model."
)
