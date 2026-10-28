"""
Plot Predicted vs. Observed pIC50 for the tuned Random Forest on the held-out Test set.

Run from the project root:
    python notebooks/19_plot_final_test_predicted_vs_observed.py

Input:
    results/tables/egfr_final_test_predictions.csv

Output:
    results/figures/egfr_tuned_test_predicted_vs_observed.png

This script only visualizes saved test predictions. It does NOT fit or tune a model.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

INPUT = Path("results/tables/egfr_final_test_predictions.csv")
OUTPUT = Path("results/figures/egfr_tuned_test_predicted_vs_observed.png")

# Candidate column names to accommodate common naming conventions.
OBSERVED_CANDIDATES = [
    "observed_pIC50"
]

PREDICTED_CANDIDATES = [
    "RandomForest_Tuned_predicted_pIC50"
]
def find_column(columns, candidates):
    lower_to_original = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    return None

if not INPUT.exists():
    raise FileNotFoundError(
        f"Could not find {INPUT}. Run notebooks/18_final_test_evaluation.py first, "
        "or update INPUT to the location of egfr_final_test_predictions.csv."
    )

df = pd.read_csv(INPUT)
observed_col = find_column(df.columns, OBSERVED_CANDIDATES)
predicted_col = find_column(df.columns, PREDICTED_CANDIDATES)

if observed_col is None or predicted_col is None:
    print("Columns found in prediction file:")
    print(df.columns.tolist())
    raise ValueError(
        "Could not identify observed and tuned-prediction columns automatically. "
        "Edit OBSERVED_CANDIDATES / PREDICTED_CANDIDATES in this script to match "
        "the column names printed above."
    )

plot_df = df[[observed_col, predicted_col]].copy()
plot_df[observed_col] = pd.to_numeric(plot_df[observed_col], errors="coerce")
plot_df[predicted_col] = pd.to_numeric(plot_df[predicted_col], errors="coerce")
plot_df = plot_df.dropna()

if plot_df.empty:
    raise ValueError("No valid numeric observed/predicted rows found in the CSV.")

observed = plot_df[observed_col].to_numpy()
predicted = plot_df[predicted_col].to_numpy()

low = min(observed.min(), predicted.min())
high = max(observed.max(), predicted.max())
padding = max((high - low) * 0.04, 0.1)
low -= padding
high += padding

fig, ax = plt.subplots(figsize=(7.5, 7))
ax.scatter(observed, predicted, alpha=0.58, s=24, edgecolors="none")
ax.plot([low, high], [low, high], linestyle="--", linewidth=1.6,
        label="Ideal prediction (y = x)")

ax.set_xlim(low, high)
ax.set_ylim(low, high)
ax.set_xlabel("Observed pIC50")
ax.set_ylabel("Predicted pIC50 (Tuned Random Forest)")
ax.set_title("Predicted vs. Observed pIC50 — Held-out Test Set")
ax.set_aspect("equal", adjustable="box")
ax.grid(True, alpha=0.25)
ax.legend(frameon=False)

n = len(plot_df)
mae = np.mean(np.abs(observed - predicted))
rmse = np.sqrt(np.mean((observed - predicted) ** 2))
ax.text(
    0.04, 0.96,
    f"n = {n}\nMAE = {mae:.3f}\nRMSE = {rmse:.3f}",
    transform=ax.transAxes, va="top", ha="left",
    bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.85, edgecolor="0.8")
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Observed column: {observed_col}")
print(f"Prediction column: {predicted_col}")
print(f"Rows plotted: {n}")
print(f"Saved figure: {OUTPUT}")
print(f"MAE from plotted rows: {mae:.4f}")
print(f"RMSE from plotted rows: {rmse:.4f}")
