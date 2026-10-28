import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load curated IC50 dataset
# -----------------------------
INPUT = "data/interim/egfr_ic50_no_missing_smiles.csv"

df = pd.read_csv(INPUT)

print("Dataset shape:", df.shape)
print("Unique compounds:", df["molecule_chembl_id"].nunique())

# Make sure activity is numeric
df["standard_value"] = pd.to_numeric(
    df["standard_value"],
    errors="coerce"
)

# Keep only valid positive IC50 values for plotting
plot_df = df[df["standard_value"] > 0].copy()

# Log10 IC50 in nM
plot_df["log10_ic50_nM"] = np.log10(plot_df["standard_value"])

print("\nIC50 summary (nM):")
print(plot_df["standard_value"].describe())

print("\nlog10(IC50) summary:")
print(plot_df["log10_ic50_nM"].describe())


# -----------------------------
# 2. Raw IC50 distribution
# -----------------------------

plt.figure(figsize=(9, 5))

# Create logarithmically spaced bins
log_bins = np.logspace(
    np.log10(plot_df["standard_value"].min()),
    np.log10(plot_df["standard_value"].max()),
    60
)

plt.hist(
    plot_df["standard_value"],
    bins=log_bins
)

plt.xlabel("IC50 (nM)")
plt.ylabel("Number of measurements")
plt.title("EGFR IC50 Distribution — Logarithmic X-axis")
plt.xscale("log")

plt.tight_layout()

plt.savefig(
    "results/figures/egfr_ic50_distribution_raw_logx.png",
    dpi=300
)

plt.show()
# -----------------------------
# 3. log10(IC50) distribution
# -----------------------------
plt.figure(figsize=(9, 5))

plt.hist(
    plot_df["log10_ic50_nM"],
    bins=60
)

plt.xlabel("log10(IC50 [nM])")
plt.ylabel("Number of measurements")
plt.title("EGFR IC50 Distribution — Log10 Scale")

plt.tight_layout()

plt.savefig(
    "results/figures/egfr_ic50_distribution_log10.png",
    dpi=300
)

plt.show()


# -----------------------------
# 4. Boxplot by assay type
# -----------------------------
assay_types = sorted(plot_df["assay_type"].dropna().unique())

data = [
    plot_df.loc[
        plot_df["assay_type"] == assay_type,
        "log10_ic50_nM"
    ]
    for assay_type in assay_types
]

plt.figure(figsize=(9, 5))

plt.boxplot(
    data,
    tick_labels=assay_types,
    showfliers=False
)

plt.xlabel("Assay Type")
plt.ylabel("log10(IC50 [nM])")
plt.title("EGFR IC50 Distribution by Assay Type")

plt.tight_layout()

plt.savefig(
    "results/figures/egfr_ic50_by_assay_type.png",
    dpi=300
)

plt.show()


# -----------------------------
# 5. Basic percentile table
# -----------------------------
percentiles = [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]

print("\nIC50 percentiles (nM):")
print(
    plot_df["standard_value"]
    .quantile(percentiles)
)

print("\nlog10(IC50) percentiles:")
print(
    plot_df["log10_ic50_nM"]
    .quantile(percentiles)
)
