import json
from collections import Counter
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. Load raw ChEMBL dataset
# ============================================================

DATA_PATH = Path("data/raw/chembl_egfr_activities.json")

with open(DATA_PATH, "r") as f:
    activities = json.load(f)

print("=" * 70)
print("INITIAL EDA — EGFR ChEMBL RAW DATASET")
print("=" * 70)


# ============================================================
# 2. Dataset overview
# ============================================================

df = pd.DataFrame(activities)

print("\n[1] DATASET OVERVIEW")
print("-" * 70)

print(f"Total activity records : {len(df):,}")

print(
    f"Unique compounds      : "
    f"{df['molecule_chembl_id'].nunique():,}"
)

print(
    f"Unique assays         : "
    f"{df['assay_chembl_id'].nunique():,}"
)

print(
    f"Unique documents      : "
    f"{df['document_chembl_id'].nunique():,}"
)


# ============================================================
# 3. Dataset structure
# ============================================================

print("\n[2] DATASET STRUCTURE")
print("-" * 70)

print(f"Number of columns: {len(df.columns)}")

print("\nColumns:")

for column in df.columns:
    print(f"- {column}")


# ============================================================
# 4. Missing-value overview
# ============================================================

print("\n[3] MISSING VALUES")
print("-" * 70)

missing = df.isna().sum()
missing_percent = (missing / len(df) * 100).round(2)

missing_table = pd.DataFrame({
    "missing_count": missing,
    "missing_percent": missing_percent
})

missing_table = missing_table[
    missing_table["missing_count"] > 0
].sort_values(
    "missing_count",
    ascending=False
)

print(missing_table)


# ============================================================
# 5. Activity type distribution
# ============================================================

print("\n[4] ACTIVITY TYPE DISTRIBUTION")
print("-" * 70)

activity_type_counts = (
    df["standard_type"]
    .value_counts(dropna=False)
)

print(activity_type_counts)

plt.figure(figsize=(10, 6))

activity_type_counts.head(15).plot(kind="bar")

plt.title("Top 15 Activity Types")
plt.xlabel("Activity Type")
plt.ylabel("Number of Records")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

plt.savefig(
    "results/figures/activity_type_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 6. IC50 subset
# ============================================================

ic50 = df[
    df["standard_type"].eq("IC50")
].copy()

print("\n[5] IC50 DATASET")
print("-" * 70)

print(f"IC50 records: {len(ic50):,}")

print(
    f"Unique IC50 compounds: "
    f"{ic50['molecule_chembl_id'].nunique():,}"
)


# ============================================================
# 7. IC50 relation distribution
# ============================================================

print("\n[6] IC50 RELATION DISTRIBUTION")
print("-" * 70)

ic50_relation = (
    ic50["standard_relation"]
    .value_counts(dropna=False)
)

print(ic50_relation)


plt.figure(figsize=(8, 5))

ic50_relation.plot(kind="bar")

plt.title("IC50 Relation Distribution")
plt.xlabel("Relation")
plt.ylabel("Number of Records")
plt.xticks(rotation=0)

plt.tight_layout()

plt.savefig(
    "results/figures/ic50_relation_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 8. IC50 unit distribution
# ============================================================

print("\n[7] IC50 UNIT DISTRIBUTION")
print("-" * 70)

ic50_units = (
    ic50["standard_units"]
    .value_counts(dropna=False)
)

print(ic50_units)


plt.figure(figsize=(8, 5))

ic50_units.plot(kind="bar")

plt.title("IC50 Unit Distribution")
plt.xlabel("Unit")
plt.ylabel("Number of Records")
plt.xticks(rotation=45, ha="right")

plt.tight_layout()

plt.savefig(
    "results/figures/ic50_unit_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 9. Clean IC50 exploratory subset
# ============================================================

clean_ic50 = ic50[
    ic50["standard_relation"].eq("=")
    & ic50["standard_units"].eq("nM")
].copy()

print("\n[8] CLEAN IC50 EXPLORATORY SUBSET")
print("-" * 70)

print(
    f"Clean IC50 measurements: "
    f"{len(clean_ic50):,}"
)

print(
    f"Unique compounds: "
    f"{clean_ic50['molecule_chembl_id'].nunique():,}"
)


# ============================================================
# 10. Measurements per compound
# ============================================================

print("\n[9] MEASUREMENTS PER COMPOUND")
print("-" * 70)

measurements_per_compound = (
    clean_ic50
    .groupby("molecule_chembl_id")
    .size()
)

print(
    f"Compounds with 1 measurement : "
    f"{(measurements_per_compound == 1).sum():,}"
)

print(
    f"Compounds with 2 measurements : "
    f"{(measurements_per_compound == 2).sum():,}"
)

print(
    f"Compounds with 3 measurements : "
    f"{(measurements_per_compound == 3).sum():,}"
)

print(
    f"Compounds with >3 measurements: "
    f"{(measurements_per_compound > 3).sum():,}"
)

print(
    f"Maximum measurements/compound: "
    f"{measurements_per_compound.max():,}"
)


plt.figure(figsize=(8, 5))

distribution = measurements_per_compound.value_counts().sort_index()

distribution.head(15).plot(kind="bar")

plt.title("Measurements per Compound")
plt.xlabel("Number of Measurements")
plt.ylabel("Number of Compounds")

plt.tight_layout()

plt.savefig(
    "results/figures/measurements_per_compound.png",
    dpi=300
)

plt.close()


# ============================================================
# 11. Assays per compound
# ============================================================

print("\n[10] ASSAYS PER COMPOUND")
print("-" * 70)

assays_per_compound = (
    clean_ic50
    .groupby("molecule_chembl_id")["assay_chembl_id"]
    .nunique()
)

print(
    f"Compounds with 1 assay : "
    f"{(assays_per_compound == 1).sum():,}"
)

print(
    f"Compounds with 2 assays: "
    f"{(assays_per_compound == 2).sum():,}"
)

print(
    f"Compounds with 3 assays: "
    f"{(assays_per_compound == 3).sum():,}"
)

print(
    f"Compounds with >3 assays: "
    f"{(assays_per_compound > 3).sum():,}"
)

print(
    f"Maximum assays/compound: "
    f"{assays_per_compound.max():,}"
)


# ============================================================
# 12. Assay type distribution
# ============================================================

print("\n[11] ASSAY TYPE DISTRIBUTION")
print("-" * 70)

assay_type = (
    clean_ic50["assay_type"]
    .value_counts(dropna=False)
)

print(assay_type)


plt.figure(figsize=(7, 5))

assay_type.plot(kind="bar")

plt.title("Assay Type Distribution")
plt.xlabel("Assay Type")
plt.ylabel("Number of IC50 Measurements")
plt.xticks(rotation=0)

plt.tight_layout()

plt.savefig(
    "results/figures/assay_type_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 13. BAO distribution
# ============================================================

print("\n[12] BAO FORMAT DISTRIBUTION")
print("-" * 70)

bao = (
    clean_ic50["bao_label"]
    .value_counts(dropna=False)
)

print(bao)


plt.figure(figsize=(9, 5))

bao.plot(kind="bar")

plt.title("BAO Format Distribution")
plt.xlabel("BAO Format")
plt.ylabel("Number of IC50 Measurements")
plt.xticks(rotation=45, ha="right")

plt.tight_layout()

plt.savefig(
    "results/figures/bao_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 14. Mutation distribution
# ============================================================

print("\n[13] MUTATION CONTEXT")
print("-" * 70)

mutation = (
    clean_ic50["assay_variant_mutation"]
    .fillna("No reported mutation")
    .value_counts()
)

print(mutation)


plt.figure(figsize=(10, 6))

mutation.head(15).plot(kind="bar")

plt.title("EGFR Assay Mutation Context")
plt.xlabel("Mutation")
plt.ylabel("Number of IC50 Measurements")
plt.xticks(rotation=45, ha="right")

plt.tight_layout()

plt.savefig(
    "results/figures/mutation_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 15. pChEMBL availability
# ============================================================

print("\n[14] pChEMBL VALUE AVAILABILITY")
print("-" * 70)

pchembl_available = clean_ic50["pchembl_value"].notna()

print(
    f"pChEMBL available: "
    f"{pchembl_available.sum():,}"
)

print(
    f"pChEMBL missing: "
    f"{(~pchembl_available).sum():,}"
)


# ============================================================
# 16. SMILES availability
# ============================================================

print("\n[15] MOLECULAR STRUCTURE AVAILABILITY")
print("-" * 70)

smiles_available = clean_ic50["canonical_smiles"].notna()

print(
    f"SMILES available: "
    f"{smiles_available.sum():,}"
)

print(
    f"SMILES missing: "
    f"{(~smiles_available).sum():,}"
)


# ============================================================
# 17. Key observations
# ============================================================

print("\n" + "=" * 70)
print("KEY OBSERVATIONS")
print("=" * 70)

print("""
1. The raw ChEMBL EGFR dataset contains multiple activity types.
   IC50 is the largest single activity type.

2. IC50 measurements contain multiple relation types and units.
   Therefore, raw IC50 records should not be used directly for modeling.

3. A preliminary exact-value IC50 subset can be defined using:
      standard_type = IC50
      standard_relation = "="
      standard_units = nM

4. Compounds can have multiple IC50 measurements and multiple assays.
   Therefore, duplicate/replicate handling will be required during curation.

5. EGFR IC50 measurements come from heterogeneous assay contexts,
   including single-protein, cell-based, and other BAO formats.

6. Mutation context is heterogeneous.
   Several EGFR variants occur in the dataset, so mutation status
   must be considered before defining the final modeling dataset.

7. Missing values are present in several fields.
   They should be handled field-by-field rather than using blanket
   row deletion.

8. This notebook performs exploratory analysis only.
   No final QSAR dataset is defined here.
""")

print("\nEDA completed.")
print("Figures saved to: results/figures/")