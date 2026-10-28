
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


# Reproducibility
RANDOM_SEED = 42

# Input
INPUT_PATH = Path("data/processed/egfr_modeling_dataset.csv")

# Outputs
OUTPUT_DIR = Path("data/processed")
QC_DIR = Path("results/tables")

SPLIT_PATH = OUTPUT_DIR / "egfr_scaffold_split.csv"
SUMMARY_PATH = QC_DIR / "egfr_scaffold_split_summary.csv"


def get_scaffold(smiles):
    """Return the Bemis-Murcko scaffold for a valid SMILES."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Invalid SMILES encountered: {smiles}")

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)

    # Acyclic molecules can have an empty Murcko scaffold.
    return Chem.MolToSmiles(scaffold, canonical=True)


def main():
    df = pd.read_csv(INPUT_PATH)

    required = {
        "molecule_chembl_id",
        "canonical_smiles",
        "pIC50",
    }
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df["molecule_chembl_id"].duplicated().any():
        raise ValueError("Duplicate molecule IDs found.")

    if df[["canonical_smiles", "pIC50"]].isna().any().any():
        raise ValueError("Missing SMILES or pIC50 values found.")

    # Assign a chemical scaffold to every molecule.
    df["scaffold"] = df["canonical_smiles"].map(get_scaffold)

    # Group by scaffold so the same scaffold cannot cross splits.
    scaffold_groups = (
        df.groupby("scaffold")["molecule_chembl_id"]
        .apply(list)
        .to_dict()
    )

    # Shuffle scaffold groups reproducibly.
    import random

    rng = random.Random(RANDOM_SEED)
    groups = list(scaffold_groups.items())
    rng.shuffle(groups)

    # Process larger scaffold groups first after shuffling.
    groups.sort(key=lambda item: len(item[1]), reverse=True)

    n_total = len(df)
    train_target = 0.80 * n_total
    valid_target = 0.10 * n_total

    train_ids = []
    valid_ids = []
    test_ids = []

    # Greedy allocation to approximate 80/10/10 by molecule count.
    for scaffold, molecule_ids in groups:
        n_group = len(molecule_ids)

        if len(train_ids) + n_group <= train_target:
            train_ids.extend(molecule_ids)
        elif len(valid_ids) + n_group <= valid_target:
            valid_ids.extend(molecule_ids)
        else:
            test_ids.extend(molecule_ids)

    split_map = {
        **{mol_id: "train" for mol_id in train_ids},
        **{mol_id: "validation" for mol_id in valid_ids},
        **{mol_id: "test" for mol_id in test_ids},
    }

    df["split"] = df["molecule_chembl_id"].map(split_map)

    if df["split"].isna().any():
        raise ValueError("Some molecules were not assigned to a split.")

    # Confirm that no scaffold appears in more than one split.
    scaffold_counts = df.groupby("scaffold")["split"].nunique()

    if (scaffold_counts > 1).any():
        raise ValueError("Scaffold leakage detected across splits.")

    # Summary for review.
    summary = (
        df.groupby("split")
        .agg(
            molecules=("molecule_chembl_id", "nunique"),
            scaffolds=("scaffold", "nunique"),
            pIC50_mean=("pIC50", "mean"),
            pIC50_std=("pIC50", "std"),
            pIC50_min=("pIC50", "min"),
            pIC50_max=("pIC50", "max"),
        )
        .reindex(["train", "validation", "test"])
        .reset_index()
    )

    summary["percentage"] = (
        summary["molecules"] / len(df) * 100
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    df.to_csv(SPLIT_PATH, index=False)
    summary.to_csv(SUMMARY_PATH, index=False)

    print("=" * 65)
    print("SCAFFOLD SPLIT SUMMARY")
    print("=" * 65)
    print(f"Total molecules: {len(df):,}")
    print(f"Unique scaffolds: {df['scaffold'].nunique():,}")
    print(f"Random seed: {RANDOM_SEED}")
    print("\nSplit summary:")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nScaffold leakage:")
    print("None" if not (scaffold_counts > 1).any() else "DETECTED")

    print("\nSaved:")
    print(SPLIT_PATH)
    print(SUMMARY_PATH)

    print("\nQC status: PASS")


if __name__ == "__main__":
    main()
