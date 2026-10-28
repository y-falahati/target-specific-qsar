from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "processed" / "egfr_pic50.csv"
OUTPUT = ROOT / "data" / "interim" / "egfr_measurement_consistency.csv"


df = pd.read_csv(INPUT)

print("=" * 70)
print("EGFR MEASUREMENT CONSISTENCY")
print("=" * 70)

print(f"\nMeasurements: {len(df):,}")
print(
    f"Compounds: "
    f"{df['molecule_chembl_id'].nunique():,}"
)


# --------------------------------------------------
# Measurements per compound
# --------------------------------------------------

counts = (
    df.groupby("molecule_chembl_id")
    .size()
    .rename("n_measurements")
)

print("\n" + "-" * 70)
print("MEASUREMENTS PER COMPOUND")
print("-" * 70)

print(counts.describe())

print("\nDistribution:")

print(
    counts.value_counts()
    .sort_index()
    .head(20)
)


# --------------------------------------------------
# Compounds with repeated measurements
# --------------------------------------------------

repeated_ids = counts[counts > 1].index

repeated = df[
    df["molecule_chembl_id"].isin(repeated_ids)
].copy()

print("\n" + "-" * 70)
print("REPEATED MEASUREMENTS")
print("-" * 70)

print(
    f"Compounds with >1 measurement: "
    f"{len(repeated_ids):,}"
)

print(
    f"Measurements involved: "
    f"{len(repeated):,}"
)


# --------------------------------------------------
# Per-compound pIC50 statistics
# --------------------------------------------------

consistency = (
    df.groupby("molecule_chembl_id")
    .agg(
        n_measurements=("pIC50", "size"),
        mean_pIC50=("pIC50", "mean"),
        median_pIC50=("pIC50", "median"),
        std_pIC50=("pIC50", "std"),
        min_pIC50=("pIC50", "min"),
        max_pIC50=("pIC50", "max"),
        n_assays=("assay_chembl_id", "nunique"),
        n_documents=("document_chembl_id", "nunique"),
    )
)

consistency["pIC50_range"] = (
    consistency["max_pIC50"]
    - consistency["min_pIC50"]
)

consistency["abs_deviation_from_median"] = (
    df.groupby("molecule_chembl_id")["pIC50"]
    .transform(
        lambda x: np.abs(x - x.median())
    )
)

# The previous transform is measurement-level,
# so calculate maximum deviation per compound.
max_dev = (
    df.assign(
        abs_dev=(
            df["pIC50"]
            - df.groupby("molecule_chembl_id")["pIC50"]
            .transform("median")
        ).abs()
    )
    .groupby("molecule_chembl_id")["abs_dev"]
    .max()
    .rename("max_abs_deviation_from_median")
)

consistency = consistency.drop(
    columns=["abs_deviation_from_median"]
)

consistency = consistency.join(max_dev)


# --------------------------------------------------
# Summary of pIC50 range
# --------------------------------------------------

print("\n" + "-" * 70)
print("pIC50 RANGE ACROSS MEASUREMENTS")
print("-" * 70)

print(
    consistency.loc[
        consistency["n_measurements"] > 1,
        "pIC50_range"
    ].describe()
)


# --------------------------------------------------
# Categories
# --------------------------------------------------

repeated_consistency = consistency[
    consistency["n_measurements"] > 1
].copy()

bins = [
    -np.inf,
    0.5,
    1.0,
    2.0,
    3.0,
    np.inf,
]

labels = [
    "<=0.5",
    "0.5-1",
    "1-2",
    "2-3",
    ">3",
]

repeated_consistency["range_category"] = pd.cut(
    repeated_consistency["pIC50_range"],
    bins=bins,
    labels=labels,
)

print("\nRange categories:")
print(
    repeated_consistency["range_category"]
    .value_counts()
    .sort_index()
)


# --------------------------------------------------
# Most heterogeneous compounds
# --------------------------------------------------

print("\n" + "-" * 70)
print("TOP COMPOUNDS BY pIC50 RANGE")
print("-" * 70)

top = (
    repeated_consistency
    .sort_values(
        "pIC50_range",
        ascending=False
    )
    .head(20)
)

print(top.to_string())


# --------------------------------------------------
# Save
# --------------------------------------------------

consistency.to_csv(
    OUTPUT
)

print("\nSaved:")
print(OUTPUT)

print("\nDone.")