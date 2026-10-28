
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


# Paths
INPUT_PATH = Path("data/processed/egfr_qsar_final.csv")
OUTPUT_PATH = Path("data/processed/egfr_rdkit_descriptors.csv")
FAILED_PATH = Path("results/tables/egfr_descriptor_failures.csv")
QC_PATH = Path("results/tables/egfr_descriptor_qc.csv")


def calculate_descriptors(smiles):
    """Calculate 2D physicochemical descriptors for one molecule."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    values = {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol),
        "RingCount": rdMolDescriptors.CalcNumRings(mol),
        "HeavyAtomCount": Lipinski.HeavyAtomCount(mol),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "AliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "HeteroatomCount": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "FormalCharge": Chem.GetFormalCharge(mol),
    }

    return values


def main():
    df = pd.read_csv(INPUT_PATH)

    descriptor_rows = []
    failures = []

    for row in df.itertuples(index=False):
        result = calculate_descriptors(row.canonical_smiles)

        if result is None:
            failures.append({
                "molecule_chembl_id": row.molecule_chembl_id,
                "canonical_smiles": row.canonical_smiles,
                "reason": "RDKit could not parse SMILES",
            })
            continue

        descriptor_rows.append({
            "molecule_chembl_id": row.molecule_chembl_id,
            **result,
        })

    descriptors = pd.DataFrame(descriptor_rows)
    failed_df = pd.DataFrame(
        failures,
        columns=["molecule_chembl_id", "canonical_smiles", "reason"],
    )

    # Save descriptors and failures for traceability.
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    FAILED_PATH.parent.mkdir(parents=True, exist_ok=True)

    descriptors.to_csv(OUTPUT_PATH, index=False)
    failed_df.to_csv(FAILED_PATH, index=False)

    feature_cols = [
        col for col in descriptors.columns
        if col != "molecule_chembl_id"
    ]

    qc = pd.DataFrame({
        "missing_n": descriptors[feature_cols].isna().sum(),
        "missing_pct": descriptors[feature_cols].isna().mean() * 100,
        "n_unique": descriptors[feature_cols].nunique(dropna=False),
        "std": descriptors[feature_cols].std(),
    })
    qc.index.name = "descriptor"
    qc.to_csv(QC_PATH)

    print("=" * 65)
    print("RDKit DESCRIPTOR GENERATION")
    print("=" * 65)
    print(f"Input compounds: {len(df):,}")
    print(f"Successful calculations: {len(descriptors):,}")
    print(f"Failed calculations: {len(failed_df):,}")
    print(f"Unique output molecule IDs: "
          f"{descriptors['molecule_chembl_id'].nunique():,}")

    print("\nDescriptor columns:")
    print(feature_cols)

    print("\nMissing values per descriptor:")
    print(descriptors[feature_cols].isna().sum().to_string())

    print("\nConstant descriptors (one unique value):")
    print(qc.index[qc["n_unique"] == 1].tolist())

    print("\nSaved:")
    print(OUTPUT_PATH)
    print(FAILED_PATH)
    print(QC_PATH)


if __name__ == "__main__":
    main()
