from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

INPUT = Path("results/tables/egfr_final_test_predictions.csv")
OUTPUT = Path("results/figures/egfr_tuned_test_residuals.png")

OBSERVED_COL = "observed_pIC50"
PREDICTED_COL = "RandomForest_Tuned_predicted_pIC50"
RESIDUAL_COL = "RandomForest_Tuned_residual"

if not INPUT.exists():
    raise FileNotFoundError(
        f"Could not find {INPUT}. Run this script from the project root."
    )

df = pd.read_csv(INPUT)
missing = [c for c in (OBSERVED_COL, PREDICTED_COL) if c not in df.columns]
if missing:
    print("Columns found:", df.columns.tolist())
    raise ValueError(f"Missing required column(s): {missing}")

observed = pd.to_numeric(df[OBSERVED_COL], errors="coerce")
predicted = pd.to_numeric(df[PREDICTED_COL], errors="coerce")
if RESIDUAL_COL in df.columns:
    residual = pd.to_numeric(df[RESIDUAL_COL], errors="coerce")
else:
    residual = observed - predicted

plot_df = pd.DataFrame({"predicted": predicted, "residual": residual})
plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()
if plot_df.empty:
    raise ValueError("No valid numeric prediction/residual rows found.")

fig, ax = plt.subplots(figsize=(8, 5.8))
ax.scatter(plot_df["predicted"], plot_df["residual"],
           alpha=0.58, s=24, edgecolors="none")
ax.axhline(0, linestyle="--", linewidth=1.6, label="Zero residual")
ax.set_xlabel("Predicted pIC50 (Tuned Random Forest)")
ax.set_ylabel("Residual (Observed − Predicted)")
ax.set_title("Residuals vs. Predicted pIC50 — Held-out Test Set")
ax.grid(True, alpha=0.25)
ax.legend(frameon=False)

mean_residual = plot_df["residual"].mean()
mae = plot_df["residual"].abs().mean()
rmse = np.sqrt(np.mean(plot_df["residual"] ** 2))
ax.text(
    0.03, 0.97,
    f"n = {len(plot_df)}\nMean residual = {mean_residual:.3f}\nMAE = {mae:.3f}\nRMSE = {rmse:.3f}",
    transform=ax.transAxes, va="top", ha="left",
    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
              alpha=0.85, edgecolor="0.8")
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Rows plotted: {len(plot_df)}")
print(f"Mean residual (observed - predicted): {mean_residual:.4f}")
print(f"MAE from plotted rows: {mae:.4f}")
print(f"RMSE from plotted rows: {rmse:.4f}")
print(f"Saved figure: {OUTPUT}")
