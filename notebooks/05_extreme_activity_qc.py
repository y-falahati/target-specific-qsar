import pandas as pd
import numpy as np

# -----------------------------
# 1. Paths
# -----------------------------
INPUT = "data/interim/egfr_ic50_equal_relation.csv"

EXTREME_OUTPUT = "data/interim/egfr_extreme_activity_qc.csv"
CLEAN_OUTPUT = "data/interim/egfr_ic50_qc_clean.csv"


# -----------------------------
# 2. Load dataset
# -----------------------------
df = pd.read_csv(INPUT)

df["standard_value"] = pd.to_numeric(
    df["standard_value"],
    errors="coerce"
)

# log10(IC50) in nM
df["log10_ic50_nM"] = np.log10(df["standard_value"])


# -----------------------------
# 3. Identify extreme records
# -----------------------------
very_low = df["log10_ic50_nM"] < -2
very_high = df["log10_ic50_nM"] > 6

extreme_mask = very_low | very_high

extreme = df[extreme_mask].copy()

extreme["extreme_group"] = np.where(
    extreme["log10_ic50_nM"] < -2,
    "very_low_IC50",
    "very_high_IC50"
)


# -----------------------------
# 4. Save extreme records
#    for traceability
# -----------------------------
extreme.to_csv(
    EXTREME_OUTPUT,
    index=False
)


# -----------------------------
# 5. Remove extreme records
# -----------------------------
clean_df = df[~extreme_mask].copy()


# -----------------------------
# 6. QC report
# -----------------------------
print("=" * 70)
print("EXTREME ACTIVITY QC — FINAL EXCLUSION")
print("=" * 70)

print(f"\nInput measurements: {len(df)}")
print(f"Input compounds: {df['molecule_chembl_id'].nunique()}")

print("\nExtreme measurements identified:")
print(f"Very low IC50 (< 0.01 nM): {very_low.sum()}")
print(f"Very high IC50 (> 1,000,000 nM): {very_high.sum()}")
print(f"Total excluded: {extreme_mask.sum()}")

print("\nAfter QC exclusion:")
print(f"Remaining measurements: {len(clean_df)}")
print(
    "Remaining compounds:",
    clean_df["molecule_chembl_id"].nunique()
)


# -----------------------------
# 7. Check that exclusion worked
# -----------------------------
remaining_low = (
    clean_df["log10_ic50_nM"] < -2
).sum()

remaining_high = (
    clean_df["log10_ic50_nM"] > 6
).sum()

print("\nRemaining extreme values:")
print(f"Very low IC50: {remaining_low}")
print(f"Very high IC50: {remaining_high}")


# -----------------------------
# 8. Check core QC fields
# -----------------------------
print("\n" + "=" * 70)
print("CORE FIELD QC")
print("=" * 70)

core_fields = [
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_value",
    "standard_units",
    "standard_relation",
    "assay_chembl_id",
]

for col in core_fields:
    print(
        f"{col:25s}: {clean_df[col].isna().sum()}"
    )


# -----------------------------
# 9. Check activity distribution
# -----------------------------
print("\n" + "=" * 70)
print("CLEAN IC50 SUMMARY")
print("=" * 70)

print(
    clean_df["standard_value"].describe()
)

print("\nClean log10(IC50) summary:")

print(
    clean_df["log10_ic50_nM"].describe()
)


# -----------------------------
# 10. Save clean dataset
# -----------------------------
clean_df.to_csv(
    CLEAN_OUTPUT,
    index=False
)

print("\n" + "=" * 70)
print(f"Saved extreme records: {EXTREME_OUTPUT}")
print(f"Saved clean dataset:   {CLEAN_OUTPUT}")
print("=" * 70)