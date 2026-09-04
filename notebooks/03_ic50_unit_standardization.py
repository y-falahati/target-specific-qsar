import json
import pandas as pd


RAW_FILE = "data/raw/chembl_egfr_activities.json"
OUTPUT_FILE = "data/interim/egfr_ic50_standardized.csv"


# --------------------------------------------------
# 1. Load raw ChEMBL data
# --------------------------------------------------

with open(RAW_FILE, "r") as f:
    data = json.load(f)

df = pd.DataFrame(data)


# --------------------------------------------------
# 2. Select IC50 measurements
# --------------------------------------------------

ic50 = df[df["standard_type"] == "IC50"].copy()

print("Total IC50 measurements:", len(ic50))


# --------------------------------------------------
# 3. Inspect original standardized units
# --------------------------------------------------

print("\nIC50 standard units:")
print(
    ic50["standard_units"]
    .value_counts(dropna=False)
)


# --------------------------------------------------
# 4. Keep only ChEMBL-standardized nM measurements
# --------------------------------------------------

ic50_nM = ic50[
    ic50["standard_units"] == "nM"
].copy()


# --------------------------------------------------
# 5. Convert standard_value to numeric
# --------------------------------------------------

ic50_nM["standard_value"] = pd.to_numeric(
    ic50_nM["standard_value"],
    errors="coerce"
)


# --------------------------------------------------
# 6. Remove records without a quantitative value
# --------------------------------------------------

before_value_filter = len(ic50_nM)

ic50_nM = ic50_nM[
    ic50_nM["standard_value"].notna()
].copy()

removed_missing_value = (
    before_value_filter - len(ic50_nM)
)


# --------------------------------------------------
# 7. Validate positive IC50 values
# --------------------------------------------------

invalid_values = (
    ic50_nM["standard_value"] <= 0
).sum()

if invalid_values > 0:
    raise ValueError(
        f"Found {invalid_values} IC50 values <= 0."
    )


# --------------------------------------------------
# 8. Final unit-standardized dataset
# --------------------------------------------------

ic50_nM.to_csv(
    OUTPUT_FILE,
    index=False
)


# --------------------------------------------------
# 9. Final report
# --------------------------------------------------

print("\n" + "=" * 50)
print("UNIT STANDARDIZATION SUMMARY")
print("=" * 50)

print(
    f"Original IC50 measurements: {len(ic50)}"
)

print(
    f"ChEMBL-standardized nM measurements: "
    f"{len(ic50_nM) + removed_missing_value}"
)

print(
    f"Removed due to missing standard_value: "
    f"{removed_missing_value}"
)

print(
    f"Final standardized IC50 measurements: "
    f"{len(ic50_nM)}"
)

print(
    f"Final unique compounds: "
    f"{ic50_nM['molecule_chembl_id'].nunique()}"
)

print(
    f"Final units: "
    f"{ic50_nM['standard_units'].unique()}"
)

print(
    f"Invalid values <= 0: "
    f"{invalid_values}"
)

print("\nSaved:")
print(OUTPUT_FILE)