
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import joblib

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


# Paths
DATA_PATH = Path("data/processed/egfr_scaffold_split.csv")
MODEL_PATH = Path("models/randomforest_baseline.joblib")
PREDICTIONS_PATH = Path(
    "results/tables/egfr_baseline_validation_predictions.csv"
)

TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

FEATURES = [
    "MolWt", "MolLogP", "TPSA", "HBD", "HBA",
    "RotatableBonds", "RingCount", "HeavyAtomCount",
    "FractionCSP3", "AromaticRings", "AliphaticRings",
    "HeteroatomCount", "FormalCharge",
]


def main():
    df = pd.read_csv(DATA_PATH)
    predictions = pd.read_csv(PREDICTIONS_PATH)

    # Select Random Forest validation predictions only.
    predictions = predictions[
        predictions["model"] == "RandomForest"
    ].copy()

    if predictions.empty:
        raise ValueError("No Random Forest validation predictions found.")

    # Merge by molecule ID to recover the chemical features.
    validation = df[df["split"] == "validation"].copy()

    analysis = predictions.merge(
        validation[["molecule_chembl_id"] + FEATURES],
        on="molecule_chembl_id",
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    if not (analysis["_merge"] == "both").all():
        raise ValueError("Some predictions could not be matched to molecules.")

    analysis = analysis.drop(columns="_merge")

    analysis["residual"] = (
        analysis["y_true_pIC50"] - analysis["y_pred_pIC50"]
    )
    analysis["absolute_error"] = analysis["residual"].abs()

    y_true = analysis["y_true_pIC50"]
    y_pred = analysis["y_pred_pIC50"]

    print("=" * 65)
    print("RANDOM FOREST VALIDATION ERROR ANALYSIS")
    print("=" * 65)
    print(f"Validation molecules: {len(analysis):,}")
    print(f"R²:   {r2_score(y_true, y_pred):.4f}")
    print(f"RMSE: {np.sqrt(mean_squared_error(y_true, y_pred)):.4f}")
    print(f"MAE:  {mean_absolute_error(y_true, y_pred):.4f}")
    print(f"Mean residual: {analysis['residual'].mean():.4f}")

    print("\nAbsolute error summary:")
    print(analysis["absolute_error"].describe().to_string())

    # Save individual predictions and errors.
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    error_path = TABLES_DIR / "egfr_randomforest_validation_errors.csv"
    analysis.to_csv(error_path, index=False)

    # Plot 1: Predicted vs observed pIC50.
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(y_true, y_pred, alpha=0.45)
    lower = min(y_true.min(), y_pred.min())
    upper = max(y_true.max(), y_pred.max())
    ax.plot([lower, upper], [lower, upper], linestyle="--")
    ax.set_xlabel("Observed pIC50")
    ax.set_ylabel("Predicted pIC50")
    ax.set_title("Random Forest: Predicted vs Observed")
    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / "egfr_rf_predicted_vs_observed.png",
        dpi=300,
    )
    plt.close(fig)

    # Plot 2: Residuals vs predicted values.
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(y_pred, analysis["residual"], alpha=0.45)
    ax.axhline(0, linestyle="--")
    ax.set_xlabel("Predicted pIC50")
    ax.set_ylabel("Residual (observed - predicted)")
    ax.set_title("Random Forest: Residual Analysis")
    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / "egfr_rf_residuals.png",
        dpi=300,
    )
    plt.close(fig)

    # Plot 3: Distribution of absolute errors.
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.hist(analysis["absolute_error"], bins=30)
    ax.set_xlabel("Absolute error in pIC50")
    ax.set_ylabel("Number of molecules")
    ax.set_title("Random Forest: Absolute Error Distribution")
    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / "egfr_rf_absolute_errors.png",
        dpi=300,
    )
    plt.close(fig)

    # Summarize errors by observed activity range.
    analysis["activity_group"] = pd.cut(
        analysis["y_true_pIC50"],
        bins=[-np.inf, 6, 7, 8, np.inf],
        labels=["<6", "6–7", "7–8", ">8"],
        right=False,
    )

    group_summary = (
        analysis.groupby("activity_group", observed=False)
        .agg(
            molecules=("molecule_chembl_id", "count"),
            MAE=("absolute_error", "mean"),
            mean_residual=("residual", "mean"),
        )
        .reset_index()
    )

    group_path = TABLES_DIR / "egfr_rf_error_by_activity.csv"
    group_summary.to_csv(group_path, index=False)

    print("\nError by observed activity:")
    print(group_summary.to_string(index=False))

    print("\nSaved:")
    print(error_path)
    print(group_path)
    print(FIGURES_DIR / "egfr_rf_predicted_vs_observed.png")
    print(FIGURES_DIR / "egfr_rf_residuals.png")
    print(FIGURES_DIR / "egfr_rf_absolute_errors.png")

    print("\nTest set was not used.")
    print("QC status: PASS")


if __name__ == "__main__":
    main()
