from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "processed" / "egfr_pic50.csv"
OUTPUT = ROOT / "data" / "processed" / "egfr_qsar_final.csv"


df = pd.read_csv(INPUT)

print("=" * 70)
print("FINAL EGFR QSAR DATASET")
print("=" * 70)

print(f"\nInput measurements: {len(df):,}")
print(
    f"Input compounds: "
    f"{df['molecule_chembl_id'].nunique():,}"
)


# --------------------------------------------------
# Aggregate measurements per compound
# --------------------------------------------------

final = (
    df.groupby("molecule_chembl_id")
    .agg(
        canonical_smiles=("canonical_smiles", "first"),

        pIC50=("pIC50", "median"),

        n_measurements=("pIC50", "size"),

        pIC50_mean=("pIC50", "mean"),
        pIC50_std=("pIC50", "std"),
        pIC50_min=("pIC50", "min"),
        pIC50_max=("pIC50", "max"),

        n_assays=("assay_chembl_id", "nunique"),
        n_documents=("document_chembl_id", "nunique"),
    )
    .reset_index()
)


# --------------------------------------------------
# Derived QC variables
# --------------------------------------------------

final["pIC50_range"] = (
    final["pIC50_max"]
    - final["pIC50_min"]
)


# --------------------------------------------------
# Validate SMILES
# --------------------------------------------------

missing_smiles = (
    final["canonical_smiles"]
    .isna()
    .sum()
)

empty_smiles = (
    final["canonical_smiles"]
    .astype(str)
    .str.strip()
    .eq("")
    .sum()
)

print("\n" + "-" * 70)
print("SMILES QC")
print("-" * 70)

print(f"Missing SMILES: {missing_smiles}")
print(f"Empty SMILES: {empty_smiles}")


# --------------------------------------------------
# Validate pIC50
# --------------------------------------------------

missing_pic50 = final["pIC50"].isna().sum()

print("\n" + "-" * 70)
print("pIC50 QC")
print("-" * 70)

print(f"Missing pIC50: {missing_pic50}")

print("\npIC50 summary:")
print(final["pIC50"].describe())


# --------------------------------------------------
# Validate uniqueness
# --------------------------------------------------

duplicate_molecules = (
    final["molecule_chembl_id"].duplicated().sum()
)

print("\n" + "-" * 70)
print("UNIQUENESS QC")
print("-" * 70)

print(
    f"Duplicate molecule IDs: "
    f"{duplicate_molecules}"
)


# --------------------------------------------------
# Distribution of measurement counts
# --------------------------------------------------

print("\n" + "-" * 70)
print("MEASUREMENTS PER COMPOUND")
print("-" * 70)

print(
    final["n_measurements"]
    .value_counts()
    .sort_index()
    .head(20)
)


# --------------------------------------------------
# High heterogeneity
# --------------------------------------------------

print("\n" + "-" * 70)
print("HIGH HETEROGENEITY")
print("-" * 70)

high_range = final[
    final["pIC50_range"] > 2
]

print(
    f"Compounds with pIC50 range > 2: "
    f"{len(high_range):,}"
)

print(
    f"Percentage: "
    f"{100 * len(high_range) / len(final):.2f}%"
)


# --------------------------------------------------
# Save
# --------------------------------------------------

final.to_csv(
    OUTPUT,
    index=False
)

print("\n" + "=" * 70)
print("FINAL DATASET SAVED")
print("=" * 70)

print(OUTPUT)

print(f"\nFinal compounds: {len(final):,}")

print("\nColumns:")
print(final.columns.tolist())

print("\nDone.")
