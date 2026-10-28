
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import joblib

from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error


# Paths
DATA_PATH = Path("data/processed/egfr_scaffold_split.csv")
MODEL_PATH = Path("models/randomforest_baseline.joblib")

TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")

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
RANDOM_SEED = 42


def main():
    df = pd.read_csv(DATA_PATH)

    validation = df[df["split"] == "validation"].copy()

    if validation.empty:
        raise ValueError("Validation dataset is empty.")

    required = FEATURES + [TARGET]
    if validation[required].isna().any().any():
        raise ValueError("Missing values found in validation data.")

    X_val = validation[FEATURES]
    y_val = validation[TARGET]

    model = joblib.load(MODEL_PATH)

    # Baseline validation performance.
    baseline_predictions = model.predict(X_val)
    baseline_rmse = mean_squared_error(
        y_val, baseline_predictions
    ) ** 0.5

    # Permute each feature repeatedly and measure the RMSE increase.
    result = permutation_importance(
        model,
        X_val,
        y_val,
        scoring="neg_root_mean_squared_error",
        n_repeats=20,
        random_state=RANDOM_SEED,
        n_jobs=-1,
    )

    importance = pd.DataFrame({
        "feature": FEATURES,
        "importance_mean": result.importances_mean,
        "importance_std": result.importances_std,
    })

    importance = importance.sort_values(
        "importance_mean", ascending=False
    ).reset_index(drop=True)

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    table_path = TABLES_DIR / "egfr_rf_permutation_importance.csv"
    figure_path = FIGURES_DIR / "egfr_rf_permutation_importance.png"

    importance.to_csv(table_path, index=False)

    # Plot from lowest to highest for a horizontal bar chart.
    plot_data = importance.sort_values("importance_mean")

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(
        plot_data["feature"],
        plot_data["importance_mean"],
        xerr=plot_data["importance_std"],
        capsize=3,
    )
    ax.axvline(0, linestyle="--")
    ax.set_xlabel("Increase in validation RMSE after permutation")
    ax.set_ylabel("Molecular descriptor")
    ax.set_title("Random Forest: Permutation Feature Importance")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

    print("=" * 65)
    print("PERMUTATION FEATURE IMPORTANCE")
    print("=" * 65)
    print(f"Validation molecules: {len(validation):,}")
    print(f"Baseline validation RMSE: {baseline_rmse:.4f}")
    print("\nFeature ranking:")
    print(importance.to_string(
        index=False,
        float_format=lambda x: f"{x:.5f}",
    ))

    print("\nSaved:")
    print(table_path)
    print(figure_path)

    print("\nTest set was not used.")
    print("QC status: PASS")


if __name__ == "__main__":
    main()
