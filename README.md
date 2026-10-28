# Target-Specific EGFR Bioactivity Prediction and QSAR Modeling

A reproducible cheminformatics and machine-learning project exploring the prediction of small-molecule bioactivity against the epidermal growth factor receptor (EGFR).

The project uses experimental bioactivity data from ChEMBL, molecular descriptors calculated with RDKit, and scikit-learn models to investigate how well molecular properties predict EGFR-associated pIC50 values.

## Project Overview

**Target:** Human epidermal growth factor receptor (EGFR)
**ChEMBL target ID:** CHEMBL203
**UniProt accession:** P00533
**Primary activity endpoint:** IC50, transformed to pIC50
**Main modeling approach:** Random Forest regression
**Data source:** ChEMBL

The goal is to build and evaluate a target-specific QSAR baseline, examine its errors, and assess its performance on compounds belonging to scaffolds excluded from training.

## Workflow

1. Retrieve EGFR-associated bioactivity records from ChEMBL.
2. Standardize IC50 units and apply data-quality filters.
3. Examine activity distributions and measurement consistency.
4. Aggregate repeated measurements at the compound level.
5. Calculate molecular descriptors using RDKit.
6. Construct the modeling dataset.
7. Split compounds into training, validation, and test sets using Bemis–Murcko scaffolds.
8. Train baseline regression models.
9. Analyze prediction errors and permutation feature importance.
10. Tune Random Forest hyperparameters using scaffold-grouped cross-validation.
11. Evaluate the selected candidate on the held-out test set.

## Dataset and Molecular Representation

The curated modeling dataset contains **6,959 compounds**.

The molecular representation uses 13 RDKit descriptors, including:

- Molecular weight (`MolWt`)
- Calculated lipophilicity (`MolLogP`)
- Topological polar surface area (`TPSA`)
- Hydrogen-bond donors and acceptors (`HBD`, `HBA`)
- Rotatable bonds and ring counts
- Heavy-atom count and heteroatom count
- Fraction of sp3-hybridized carbon atoms (`FractionCSP3`)
- Formal charge

The endpoint is defined as:

`pIC50 = -log10(IC50 in molar units)`

A difference of one pIC50 unit corresponds to a tenfold difference in IC50 concentration.

Repeated measurements were aggregated at the compound level using the median pIC50. The curation workflow applied defined unit, relation, activity-range, assay-format, and mutation-reporting filters. The absence of a reported mutation should not be interpreted as independent proof that every record represents a confirmed wild-type assay.

## Data Splitting and Evaluation

A Bemis–Murcko scaffold split was used to reduce structural leakage between the training, validation, and test sets.

| Split | Compounds | Approximate proportion |
|---|---:|---:|
| Training | 5,567 | 80% |
| Validation | 695 | 10% |
| Test | 697 | 10% |
| **Total** | **6,959** | **100%** |

The split was generated with random seed 42. No Bemis–Murcko scaffold was shared across the three partitions.

### Held-out Test Results

The following metrics were obtained from the final test set, which was not used for model fitting or hyperparameter selection.

| Model | R² | RMSE | MAE |
|---|---:|---:|---:|
| Random Forest baseline | 0.251 | 1.131 | 0.863 |
| Tuned Random Forest | 0.285 | 1.105 | 0.860 |

The tuned model improved R² and RMSE modestly compared with the baseline. The improvement in MAE was small, so the results should not be interpreted as evidence of high predictive accuracy.

The test results represent performance on this particular scaffold split. External validation on an independent dataset has not yet been established.

## Key Findings

- A Random Forest model performed better than the linear regression baseline on the validation set.
- Scaffold-based splitting provides a more demanding evaluation than a random compound split when structural generalization is important.
- Hyperparameter tuning produced a modest improvement in test-set R² and RMSE.
- Residual analysis indicated a tendency toward average overprediction of pIC50.
- Permutation importance identified several descriptors associated with model performance, but importance does not establish a causal relationship between a descriptor and biological activity.

## Repository Structure

```text
target-specific-qsar/
├── data/
│   ├── raw/
│   ├── interim/
│   └── processed/
├── docs/
│   └── data_dictionary.md
├── models/
├── notebooks/
│   ├── 01_chembl_data_acquisition.py
│   ├── 02_initial_eda.py
│   ├── 03_ic50_unit_standardization.py
│   └── ...
├── reports/
│   ├── egfr_final_test_evaluation_report.md
│   └── egfr_model_limitations.md
├── results/
│   ├── figures/
│   └── tables/
├── environment.yml
├── requirements.txt
├── .gitignore
└── README.md


Raw and processed datasets, generated result tables and figures, and serialized model files are excluded from Git by default. Run the workflow locally to regenerate the relevant artifacts.
Installation
This project uses Conda and Python 3.11.
git clone https://github.com/y-falahati/target-specific-qsar.git
cd target-specific-qsar

conda env create -f environment.yml
conda activate molecular-activity

python -c "import rdkit, numpy, pandas, sklearn, matplotlib, seaborn, joblib; print('Core imports: OK')"


Running the Workflow
The numbered Python scripts in notebooks/ document the analysis stages. Run them in sequence after reviewing their input paths and prerequisites.

For example:
python notebooks/14_train_baseline_models.py
python notebooks/15_analyze_baseline_errors.py
python notebooks/16_feature_importance.py
python notebooks/17_tune_random_forest.py
python notebooks/18_final_test_evaluation.py
python notebooks/19_plot_final_test_predicted_vs_observed.py
python notebooks/20_plot_final_test_residuals.py


The later scripts depend on datasets and model files generated by earlier stages. The input bioactivity data must be acquired or prepared before running the full pipeline.
Limitations
- The model uses a compact set of 13 molecular descriptors rather than a comprehensive representation of molecular structure.
- Similar descriptor values do not guarantee similar biological activity.
- Performance is moderate, and prediction errors remain substantial for some compounds.
- Experimental assay heterogeneity and compound-level aggregation can introduce uncertainty.
- Permutation importance describes model reliance, not biological causality.
- Applicability-domain assessment and independent external validation remain important next steps.
- The current model is intended for computational exploration and compound prioritization, not as a substitute for experimental testing.
Future Directions
- Compare descriptor-based models with molecular fingerprints.
- Evaluate additional scaffold-aware modeling strategies.
- Characterize the applicability domain and prediction uncertainty.
- Validate performance on an independent, carefully curated dataset.
- Extend the workflow toward virtual screening after validation.
Acknowledgments
Bioactivity records are sourced from ChEMBL. Molecular parsing and descriptor calculation use RDKit, and model development uses scikit-learn.
Disclaimer
This repository is a research and learning project. Predictions are computational estimates and should not be interpreted as experimentally confirmed activity, clinical evidence, or therapeutic recommendations.