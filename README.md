# Target-Specific QSAR & Virtual Screening

## Overview

This project develops a target-specific QSAR pipeline for predicting molecular bioactivity and prioritizing compounds for virtual screening.

The project uses experimental bioactivity data from ChEMBL and molecular representations generated with RDKit.

## Objective

The main objective is to investigate whether molecular structure can be used to predict the bioactivity of compounds against a specific therapeutic target and to use the resulting model for computational compound prioritization.

## Project Workflow

ChEMBL
→ Target Selection
→ Bioactivity Acquisition
→ Data Curation
→ Molecular Standardization
→ Molecular Descriptors & Fingerprints
→ QSAR Modeling
→ Model Evaluation
→ Scaffold Split Evaluation
→ Applicability Domain
→ Virtual Screening
→ Compound Prioritization

## Project Structure

```text
target-specific-qsar/
├── data/
│   ├── raw/
│   ├── interim/
│   └── processed/
├── notebooks/
├── src/
├── results/
│   ├── figures/
│   └── tables/
├── README.md
├── environment.yml
├── requirements.txt
└── .gitignore
```
## Technologies:

Python

RDKit

pandas

NumPy

scikit-learn

hEMBL

Git / GitHub

----

