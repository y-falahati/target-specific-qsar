import json
from pathlib import Path

import pandas as pd


# -----------------------------
# Paths
# -----------------------------
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "interim" / "egfr_ic50_qc_clean.csv"
OUT_DIR = ROOT / "data" / "interim"

df = pd.read_csv(INPUT)

print("=" * 70)
print("EGFR BIOLOGICAL SCOPE PROFILING")
print("=" * 70)

print(f"\nInput measurements: {len(df)}")
print(f"Input compounds: {df['molecule_chembl_id'].nunique()}")
print(f"Input assays: {df['assay_chembl_id'].nunique()}")


# -----------------------------
# 1. Assay type
# -----------------------------
print("\n" + "-" * 70)
print("1. ASSAY TYPE")
print("-" * 70)

assay_type_summary = (
    df.groupby("assay_type")
    .agg(
        measurements=("activity_id", "count"),
        compounds=("molecule_chembl_id", "nunique"),
        assays=("assay_chembl_id", "nunique"),
    )
    .sort_values("measurements", ascending=False)
)

print(assay_type_summary)


# -----------------------------
# 2. BAO label
# -----------------------------
print("\n" + "-" * 70)
print("2. BAO LABEL")
print("-" * 70)

bao_summary = (
    df.groupby("bao_label")
    .agg(
        measurements=("activity_id", "count"),
        compounds=("molecule_chembl_id", "nunique"),
        assays=("assay_chembl_id", "nunique"),
    )
    .sort_values("measurements", ascending=False)
)

print(bao_summary)


# -----------------------------
# 3. Mutation context
# -----------------------------
print("\n" + "-" * 70)
print("3. MUTATION CONTEXT")
print("-" * 70)

mutation_summary = (
    df["assay_variant_mutation"]
    .fillna("NOT_REPORTED")
    .value_counts()
)

print(mutation_summary.head(30))


# -----------------------------
# 4. Cross-tab: assay type x BAO
# -----------------------------
print("\n" + "-" * 70)
print("4. ASSAY TYPE x BAO LABEL")
print("-" * 70)

cross_bao = pd.crosstab(
    df["assay_type"],
    df["bao_label"],
)

print(cross_bao)


# -----------------------------
# 5. Candidate biochemical scope
# -----------------------------
print("\n" + "-" * 70)
print("5. CANDIDATE BIOCHEMICAL SCOPE")
print("-" * 70)

candidate = df[
    (df["bao_label"] == "single protein format")
    & (
        df["assay_variant_mutation"].isna()
        | (df["assay_variant_mutation"].str.strip() == "")
    )
].copy()

print(f"Candidate measurements: {len(candidate)}")
print(
    f"Candidate compounds: "
    f"{candidate['molecule_chembl_id'].nunique()}"
)
print(
    f"Candidate assays: "
    f"{candidate['assay_chembl_id'].nunique()}"
)


# -----------------------------
# 6. Inspect assay descriptions
# -----------------------------
print("\n" + "-" * 70)
print("6. CANDIDATE ASSAY DESCRIPTIONS")
print("-" * 70)

description_summary = (
    candidate[
        [
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "bao_label",
            "assay_variant_mutation",
        ]
    ]
    .drop_duplicates()
    .sort_values("assay_chembl_id")
)

print(
    description_summary.head(100).to_string(index=False)
)


# -----------------------------
# 7. Search descriptions for WT / mutant language
# -----------------------------
print("\n" + "-" * 70)
print("7. ASSAY DESCRIPTION KEYWORDS")
print("-" * 70)

description = (
    df["assay_description"]
    .fillna("")
    .astype(str)
    .str.lower()
)

keywords = [
    "wild type",
    "wild-type",
    "wt egfr",
    "mutant",
    "mutation",
    "l858r",
    "t790m",
    "c797s",
    "l861q",
]

for keyword in keywords:
    mask = description.str.contains(keyword, regex=False)

    print(
        f"{keyword:15s}: "
        f"{mask.sum():5d} measurements | "
        f"{df.loc[mask, 'molecule_chembl_id'].nunique():5d} compounds | "
        f"{df.loc[mask, 'assay_chembl_id'].nunique():5d} assays"
    )


# -----------------------------
# 8. Candidate descriptions containing WT language
# -----------------------------
print("\n" + "-" * 70)
print("8. WT-LIKE CANDIDATE ASSAYS")
print("-" * 70)

wt_mask = description.str.contains(
    "wild type|wild-type|wt egfr",
    regex=True,
)

wt_candidate = df[
    wt_mask
    & (df["bao_label"] == "single protein format")
].copy()

wt_assays = (
    wt_candidate[
        [
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "bao_label",
            "assay_variant_mutation",
        ]
    ]
    .drop_duplicates()
    .sort_values("assay_chembl_id")
)

print(f"WT-like measurements: {len(wt_candidate)}")
print(
    f"WT-like compounds: "
    f"{wt_candidate['molecule_chembl_id'].nunique()}"
)
print(
    f"WT-like assays: "
    f"{wt_candidate['assay_chembl_id'].nunique()}"
)

print("\nExamples:")
print(wt_assays.head(100).to_string(index=False))


# -----------------------------
# 9. Save summaries
# -----------------------------
assay_type_summary.to_csv(
    OUT_DIR / "egfr_scope_assay_type_summary.csv"
)

bao_summary.to_csv(
    OUT_DIR / "egfr_scope_bao_summary.csv"
)

mutation_summary.rename("measurements").to_csv(
    OUT_DIR / "egfr_scope_mutation_summary.csv"
)

cross_bao.to_csv(
    OUT_DIR / "egfr_scope_assay_type_x_bao.csv"
)

wt_assays.to_csv(
    OUT_DIR / "egfr_scope_wt_like_assays.csv",
    index=False,
)

print("\n" + "=" * 70)
print("Scope profiling complete.")
print("=" * 70)