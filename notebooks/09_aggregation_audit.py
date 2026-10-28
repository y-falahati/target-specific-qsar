from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "processed" / "egfr_pic50.csv"


df = pd.read_csv(INPUT)

print("=" * 70)
print("EGFR AGGREGATION AUDIT")
print("=" * 70)


# --------------------------------------------------
# Compound-level structure
# --------------------------------------------------

compound = (
    df.groupby("molecule_chembl_id")
    .agg(
        n_measurements=("pIC50", "size"),
        n_assays=("assay_chembl_id", "nunique"),
        n_documents=("document_chembl_id", "nunique"),
        n_assay_types=("assay_type", "nunique"),
        median_pIC50=("pIC50", "median"),
        min_pIC50=("pIC50", "min"),
        max_pIC50=("pIC50", "max"),
    )
)

compound["pIC50_range"] = (
    compound["max_pIC50"]
    - compound["min_pIC50"]
)


# --------------------------------------------------
# Repeated compounds
# --------------------------------------------------

repeated = compound[
    compound["n_measurements"] > 1
].copy()

print("\n" + "-" * 70)
print("REPEATED COMPOUNDS")
print("-" * 70)

print(f"Repeated compounds: {len(repeated):,}")


# --------------------------------------------------
# Assay diversity
# --------------------------------------------------

print("\n" + "-" * 70)
print("ASSAY DIVERSITY")
print("-" * 70)

print(
    repeated["n_assays"]
    .describe()
)

print(
    "\nCompounds with >1 assay:",
    (repeated["n_assays"] > 1).sum()
)

print(
    "Compounds with exactly 1 assay:",
    (repeated["n_assays"] == 1).sum()
)


# --------------------------------------------------
# Document diversity
# --------------------------------------------------

print("\n" + "-" * 70)
print("DOCUMENT DIVERSITY")
print("-" * 70)

print(
    repeated["n_documents"]
    .describe()
)

print(
    "\nCompounds with >1 document:",
    (repeated["n_documents"] > 1).sum()
)


# --------------------------------------------------
# Assay type diversity
# --------------------------------------------------

print("\n" + "-" * 70)
print("ASSAY TYPE DIVERSITY")
print("-" * 70)

print(
    repeated["n_assay_types"]
    .value_counts()
    .sort_index()
)

print(
    "\nCompounds with >1 assay type:",
    (repeated["n_assay_types"] > 1).sum()
)


# --------------------------------------------------
# pIC50 range
# --------------------------------------------------

print("\n" + "-" * 70)
print("pIC50 RANGE")
print("-" * 70)

print(
    repeated["pIC50_range"]
    .describe()
)


# --------------------------------------------------
# Important cases:
# same assay vs multiple assays
# --------------------------------------------------

same_assay = repeated[
    repeated["n_assays"] == 1
]

multiple_assays = repeated[
    repeated["n_assays"] > 1
]

print("\n" + "-" * 70)
print("pIC50 RANGE BY ASSAY STRUCTURE")
print("-" * 70)

print("\nExactly one assay:")
print(
    same_assay["pIC50_range"].describe()
)

print("\nMultiple assays:")
print(
    multiple_assays["pIC50_range"].describe()
)


# --------------------------------------------------
# Potentially problematic compounds
# --------------------------------------------------

print("\n" + "-" * 70)
print("HIGH HETEROGENEITY")
print("-" * 70)

high_heterogeneity = repeated[
    repeated["pIC50_range"] > 2
].sort_values(
    "pIC50_range",
    ascending=False
)

print(
    f"Compounds with pIC50 range > 2: "
    f"{len(high_heterogeneity):,}"
)

print(
    high_heterogeneity.head(30).to_string()
)


# --------------------------------------------------
# Save
# --------------------------------------------------

output = (
    ROOT
    / "data"
    / "interim"
    / "egfr_aggregation_audit.csv"
)

compound.to_csv(output)

print("\nSaved:")
print(output)

print("\nDone.")