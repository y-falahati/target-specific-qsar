from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------
# Paths
# --------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "interim" / "egfr_ic50_qc_clean.csv"

OUTPUT_SCOPE = (
    ROOT / "data" / "processed" / "egfr_ic50_scope.csv"
)

OUTPUT_PIC50 = (
    ROOT / "data" / "processed" / "egfr_pic50.csv"
)


# --------------------------------------------------
# Load
# --------------------------------------------------
df = pd.read_csv(INPUT)

print("=" * 70)
print("EGFR BIOLOGICAL SCOPE + pIC50")
print("=" * 70)

print(f"\nInput measurements: {len(df):,}")
print(
    f"Input compounds: "
    f"{df['molecule_chembl_id'].nunique():,}"
)


# --------------------------------------------------
# Biological scope
# --------------------------------------------------
# Scope:
#   - single protein format
#   - no reported mutation
#   - IC50
#   - exact relation (=)
#
# IC50 and relation were already enforced upstream,
# but we keep the checks explicit here.
# --------------------------------------------------

scope_mask = (
    (df["standard_type"] == "IC50")
    & (df["standard_relation"] == "=")
    & (df["standard_units"] == "nM")
    & (df["bao_label"] == "single protein format")
    & (
        df["assay_variant_mutation"].isna()
        | (
            df["assay_variant_mutation"]
            .astype(str)
            .str.strip()
            .eq("")
        )
    )
)

scope = df.loc[scope_mask].copy()


# --------------------------------------------------
# Scope summary
# --------------------------------------------------
print("\n" + "-" * 70)
print("BIOLOGICAL SCOPE")
print("-" * 70)

print(f"Scope measurements: {len(scope):,}")
print(
    f"Scope compounds: "
    f"{scope['molecule_chembl_id'].nunique():,}"
)
print(
    f"Scope assays: "
    f"{scope['assay_chembl_id'].nunique():,}"
)
print(
    f"Scope documents: "
    f"{scope['document_chembl_id'].nunique():,}"
)

print("\nAssay types:")
print(scope["assay_type"].value_counts())

print("\nBAO labels:")
print(scope["bao_label"].value_counts())

print("\nMutation field:")
print(
    scope["assay_variant_mutation"]
    .fillna("NOT_REPORTED")
    .value_counts()
)


# --------------------------------------------------
# Validate IC50 values
# --------------------------------------------------
scope["standard_value"] = pd.to_numeric(
    scope["standard_value"],
    errors="coerce",
)

missing_values = scope["standard_value"].isna().sum()
nonpositive_values = (
    scope["standard_value"] <= 0
).sum()

print("\n" + "-" * 70)
print("IC50 VALIDATION")
print("-" * 70)

print(f"Missing IC50 values: {missing_values}")
print(f"Non-positive IC50 values: {nonpositive_values}")


if missing_values > 0:
    raise ValueError(
        "Missing standard_value found in scope."
    )

if nonpositive_values > 0:
    raise ValueError(
        "Non-positive IC50 found in scope."
    )


# --------------------------------------------------
# Convert nM -> M
# --------------------------------------------------
scope["ic50_M"] = (
    scope["standard_value"] * 1e-9
)


# --------------------------------------------------
# Calculate pIC50
# --------------------------------------------------
scope["pIC50"] = -np.log10(
    scope["ic50_M"]
)


# --------------------------------------------------
# Numerical validation
# --------------------------------------------------
if not np.isfinite(scope["pIC50"]).all():
    raise ValueError(
        "Non-finite pIC50 values detected."
    )


# --------------------------------------------------
# pIC50 summary
# --------------------------------------------------
print("\n" + "-" * 70)
print("pIC50 SUMMARY")
print("-" * 70)

print(
    scope["pIC50"]
    .describe()
)

print("\nIC50 [nM] percentiles:")
print(
    scope["standard_value"]
    .quantile([0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
)

print("\npIC50 percentiles:")
print(
    scope["pIC50"]
    .quantile([0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
)


# --------------------------------------------------
# Check mathematical consistency
# --------------------------------------------------
check_pic50 = -np.log10(
    scope["standard_value"] * 1e-9
)

max_difference = np.abs(
    scope["pIC50"] - check_pic50
).max()

print("\nMaximum pIC50 calculation difference:")
print(max_difference)

if max_difference > 1e-12:
    raise ValueError(
        "pIC50 calculation consistency check failed."
    )


# --------------------------------------------------
# Save scope dataset
# --------------------------------------------------
scope.to_csv(
    OUTPUT_SCOPE,
    index=False,
)

# A compact ML-ready table
pic50_columns = [
    "activity_id",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "canonical_smiles",
    "assay_chembl_id",
    "document_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_units",
    "standard_value",
    "ic50_M",
    "pIC50",
    "assay_type",
    "bao_label",
    "assay_variant_mutation",
    "potential_duplicate",
]

pic50 = scope[pic50_columns].copy()

pic50.to_csv(
    OUTPUT_PIC50,
    index=False,
)


# --------------------------------------------------
# Final report
# --------------------------------------------------
print("\n" + "=" * 70)
print("FINAL OUTPUT")
print("=" * 70)

print(f"Scope dataset:")
print(OUTPUT_SCOPE)

print(f"\npIC50 dataset:")
print(OUTPUT_PIC50)

print("\nFinal measurements:")
print(f"{len(pic50):,}")

print("Final compounds:")
print(
    f"{pic50['molecule_chembl_id'].nunique():,}"
)

print("\nColumns:")
print(
    pic50.columns.tolist()
)

print("\nDone.")