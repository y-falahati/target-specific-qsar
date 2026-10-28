
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------
# Configuration
# --------------------------------------------------
RANDOM_SEED = 42

INPUT_PATH = Path("data/processed/egfr_scaffold_split.csv")
RESULTS_DIR = Path("results/tables")
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
    """Calculate regression evaluation metrics."""
    return {
        "R2": r2_score(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "MAE": mean_absolute_error(y_true, y_pred),
    }


def main():
    df = pd.read_csv(INPUT_PATH)

    required = FEATURES + [TARGET, "split", "molecule_chembl_id"]
    missing = set(required) - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df[required].isna().any().any():
        raise ValueError("Missing values detected in modeling columns.")

    # Separate the predefined scaffold splits.
    train = df[df["split"] == "train"].copy()
    validation = df[df["split"] == "validation"].copy()

    if train.empty or validation.empty:
        raise ValueError("Train or validation split is empty.")

    X_train = train[FEATURES]
    y_train = train[TARGET]

    X_val = validation[FEATURES]
    y_val = validation[TARGET]

    # Define baseline models.
    models = {
        "LinearRegression": Pipeline([
            ("scaler", StandardScaler()),
            ("model", LinearRegression()),
        ]),
        "RandomForest": RandomForestRegressor(
            n_estimators=300,
            min_samples_leaf=2,
            random_state=RANDOM_SEED,
            n_jobs=-1,
        ),
    }

    results = []
    prediction_frames = []

    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for model_name, model in models.items():
        print(f"\nTraining: {model_name}")

        # Fit using TRAIN only.
        model.fit(X_train, y_train)

        # Compare models using VALIDATION only.
        val_predictions = model.predict(X_val)

        metrics = calculate_metrics(y_val, val_predictions)

        results.append({
            "model": model_name,
            "train_n": len(train),
            "validation_n": len(validation),
            **metrics,
        })

        prediction_frames.append(pd.DataFrame({
            "molecule_chembl_id": validation["molecule_chembl_id"].values,
            "y_true_pIC50": y_val.values,
            "y_pred_pIC50": val_predictions,
            "model": model_name,
        }))

        # Save fitted model for reproducibility.
        model_path = MODELS_DIR / f"{model_name.lower()}_baseline.joblib"
        joblib.dump(model, model_path)

        print(f"R²:   {metrics['R2']:.4f}")
        print(f"RMSE: {metrics['RMSE']:.4f}")
        print(f"MAE:  {metrics['MAE']:.4f}")
        print(f"Saved model: {model_path}")

    results_df = pd.DataFrame(results).sort_values(
        "RMSE", ascending=True
    )

    predictions_df = pd.concat(
        prediction_frames, ignore_index=True
    )

    results_path = RESULTS_DIR / "egfr_baseline_validation_metrics.csv"
    predictions_path = (
        RESULTS_DIR / "egfr_baseline_validation_predictions.csv"
    )

    results_df.to_csv(results_path, index=False)
    predictions_df.to_csv(predictions_path, index=False)

    print("\n" + "=" * 65)
    print("BASELINE MODEL COMPARISON — VALIDATION")
    print("=" * 65)
    print(results_df.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    ))

    print("\nSaved:")
    print(results_path)
    print(predictions_path)

    print("\nTest set was not used.")
    print("QC status: PASS")


if __name__ == "__main__":
    main()
