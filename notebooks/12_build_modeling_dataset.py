
from pathlib import Path

import pandas as pd


# Input files
ACTIVITY_PATH = Path("data/processed/egfr_qsar_final.csv")
DESCRIPTORS_PATH = Path("data/processed/egfr_rdkit_descriptors.csv")

# Output files
OUTPUT_PATH = Path("data/processed/egfr_modeling_dataset.csv")
QC_PATH = Path("results/tables/egfr_modeling_dataset_qc.csv")


def main():
    activity = pd.read_csv(ACTIVITY_PATH)
    descriptors = pd.read_csv(DESCRIPTORS_PATH)

    # Check required columns
    required_activity = {
        "molecule_chembl_id",
        "canonical_smiles",
        "pIC50",
    }
    required_descriptors = {"molecule_chembl_id"}

    missing_activity = required_activity - set(activity.columns)
    missing_descriptors = required_descriptors - set(descriptors.columns)

    if missing_activity:
        raise ValueError(
            f"Missing activity columns: {sorted(missing_activity)}"
        )

    if missing_descriptors:
        raise ValueError(
            f"Missing descriptor columns: {sorted(missing_descriptors)}"
        )

    # Check molecule IDs before merging
    if activity["molecule_chembl_id"].duplicated().any():
        raise ValueError("Duplicate molecule IDs in activity dataset.")

    if descriptors["molecule_chembl_id"].duplicated().any():
        raise ValueError("Duplicate molecule IDs in descriptor dataset.")

    activity_ids = set(activity["molecule_chembl_id"])
    descriptor_ids = set(descriptors["molecule_chembl_id"])

    only_activity = activity_ids - descriptor_ids
    only_descriptors = descriptor_ids - activity_ids

    # Merge using molecule ID, preserving activity rows
    merged = activity.merge(
        descriptors,
        on="molecule_chembl_id",
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    unmatched = merged.loc[
        merged["_merge"] != "both",
        ["molecule_chembl_id", "_merge"],
    ]

    merged = merged.drop(columns="_merge")

    # Identify descriptor columns
    descriptor_cols = [
        col for col in descriptors.columns
        if col != "molecule_chembl_id"
    ]

    # Dataset quality checks
    missing_descriptor_cells = int(
        merged[descriptor_cols].isna().sum().sum()
    )

    qc_rows = [
        ("activity_rows", len(activity)),
        ("descriptor_rows", len(descriptors)),
        ("merged_rows", len(merged)),
        ("unique_molecule_ids", merged["molecule_chembl_id"].nunique()),
        ("activity_ids_without_descriptors", len(only_activity)),
        ("descriptor_ids_without_activity", len(only_descriptors)),
        ("unmatched_activity_rows", len(unmatched)),
        ("missing_pIC50", int(merged["pIC50"].isna().sum())),
        ("missing_SMILES", int(merged["canonical_smiles"].isna().sum())),
        ("missing_descriptor_cells", missing_descriptor_cells),
        ("descriptor_count", len(descriptor_cols)),
    ]

    qc = pd.DataFrame(qc_rows, columns=["check", "value"])

    # Save outputs
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(OUTPUT_PATH, index=False)
    qc.to_csv(QC_PATH, index=False)

    print("=" * 65)
    print("QSAR MODELING DATASET")
    print("=" * 65)
    print(f"Activity rows: {len(activity):,}")
    print(f"Descriptor rows: {len(descriptors):,}")
    print(f"Merged rows: {len(merged):,}")
    print(f"Unique molecules: {merged['molecule_chembl_id'].nunique():,}")
    print(f"Descriptor features: {len(descriptor_cols)}")
    print(f"Activity IDs without descriptors: {len(only_activity)}")
    print(f"Descriptor IDs without activity: {len(only_descriptors)}")
    print(f"Unmatched activity rows: {len(unmatched)}")
    print(f"Missing pIC50 values: {merged['pIC50'].isna().sum()}")
    print(f"Missing SMILES values: {merged['canonical_smiles'].isna().sum()}")
    print(f"Missing descriptor cells: {missing_descriptor_cells}")

    print("\nTarget summary:")
    print(merged["pIC50"].describe().to_string())

    print("\nSaved:")
    print(OUTPUT_PATH)
    print(QC_PATH)

    if (
        len(merged) != len(activity)
        or merged["molecule_chembl_id"].nunique() != len(activity)
        or len(unmatched) > 0
        or missing_descriptor_cells > 0
        or merged["pIC50"].isna().any()
        or merged["canonical_smiles"].isna().any()
    ):
        raise ValueError(
            "Modeling dataset QC failed. Review the checks above."
        )

    print("\nQC status: PASS")


if __name__ == "__main__":
    main()
