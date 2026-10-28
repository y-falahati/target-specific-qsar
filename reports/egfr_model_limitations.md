import pypandoc
from pathlib import Path

content = r"""
# Model Limitations

## 1. Limited molecular representation

The tuned Random Forest model uses 13 RDKit-computed molecular descriptors: `MolWt`, `MolLogP`, `TPSA`, `HBD`, `HBA`, `RotatableBonds`, `RingCount`, `HeavyAtomCount`, `FractionCSP3`, `AromaticRings`, `AliphaticRings`, `HeteroatomCount`, and `FormalCharge`.

These descriptors summarize selected physicochemical properties and broad aspects of molecular composition. They do not encode the full molecular graph or all local chemical environments. Consequently, compounds with different atom connectivity, substitution patterns, stereochemistry, or functional-group arrangements may receive similar descriptor values even when their EGFR activity differs. Conversely, similar overall descriptor profiles do not guarantee similar binding modes or biological activity.

This is a limitation of the current feature representation, not proof that the model cannot learn any structure–activity relationships.

## 2. Structural diversity and scaffold generalization

The dataset was split using Bemis–Murcko scaffolds to reduce scaffold leakage between training, validation, and test sets. This is a more demanding evaluation than a random compound-level split, but it does not guarantee reliable predictions for every novel scaffold or chemical series.

The model may perform unevenly across scaffold families, particularly when a test compound belongs to a chemical region that is sparsely represented in the training data. Overall test metrics can conceal such scaffold-specific errors.

## 3. Moderate held-out test performance

On the held-out test set (697 compounds), the tuned Random Forest achieved:

| Metric | Test result |
|---|---:|
| R² | 0.2847 |
| RMSE | 1.1053 |
| MAE | 0.8595 |
| Mean residual (observed − predicted) | −0.1904 |

These results indicate useful but limited predictive signal, rather than high-precision activity prediction. The negative mean residual indicates a modest average tendency to overpredict pIC50 under the project's residual definition. Individual predictions can deviate substantially from measured values.

Because pIC50 is a base-10 logarithmic transformation of molar potency, an error of 1 pIC50 unit corresponds to a tenfold difference in concentration. Therefore, errors around one pIC50 unit can be biologically meaningful.

## 4. Uneven errors and outliers

The test residual plot shows substantial spread around the zero-residual line, including some large positive and negative errors. A mean residual close to zero would not, by itself, establish accurate predictions because positive and negative errors can cancel each other out. MAE, RMSE, residual plots, and examination of individual high-error compounds should be considered together.

The residual plot is diagnostic, not evidence that any particular cause of error has been established. Potential explanations—such as underrepresented scaffolds, specific interactions not captured by the descriptors, assay variability, or limitations in the curated activity data—require further investigation.

## 5. Activity data and aggregation

The model target is based on curated ChEMBL IC50 measurements transformed to pIC50 and aggregated at the compound level using the median. Aggregation gives each compound a single target value, but can hide disagreement among measurements, assays, documents, and experimental conditions. IC50 values may also depend on assay design and conditions; they should not be assumed to be perfectly comparable across all records.

The curation filters and the absence of a reported mutation do not prove that every assay represents experimentally confirmed wild-type EGFR under identical conditions.

## 6. Feature importance is not causality

Permutation importance describes how much the fitted model's validation performance depends on a feature under the evaluated setup. It does not establish that a descriptor causes EGFR inhibition. Correlated descriptors can share or mask importance. The zero permutation importance observed for `FormalCharge` means that shuffling it did not measurably reduce validation performance in this setup; it does not prove formal charge is biologically irrelevant.

## 7. Applicability domain and uncertainty

The current workflow does not provide a calibrated uncertainty estimate or a validated applicability-domain decision for each candidate. A prediction for a molecule unlike the training compounds should therefore be treated cautiously. The model should not be assumed to be equally reliable across all chemical structures or activity ranges.

## 8. Test-set integrity and scope of conclusions

The held-out test set was used for final evaluation after model selection and tuning on training/validation data. Its results should be reported as the final evaluation for this split. Repeatedly using the test set to select features, tune hyperparameters, or choose among new models would make it part of the model-development process and weaken its role as an independent holdout.

The test set is still drawn from the same curated project dataset; it is not an external validation study. Performance on an independent dataset, with carefully harmonized target and assay definitions, remains to be established.

## 9. Practical interpretation

This model is best treated as a baseline virtual-screening or prioritization tool, not as a replacement for experimental testing. Its predictions can help rank compounds for follow-up, but should be combined with chemical review, scaffold coverage, data-quality checks, and experimental validation. Predictions for unfamiliar scaffolds or compounds with large expected errors warrant particular caution.

## Summary

The main limitations are the compact 13-descriptor representation, potentially different biological behavior among structurally distinct compounds with similar descriptor profiles, moderate held-out test performance, heterogeneous experimental activity data, uneven residuals, and the absence of validated per-compound uncertainty and external validation. These limitations define the current model's scope and should guide future improvements, such as richer molecular representations, scaffold-aware error analysis, uncertainty estimation, and external validation—without reusing the current test set for further tuning.
"""
out = Path("/mnt/data/egfr_model_limitations.md")
pypandoc.convert_text(content, "md", format="md", outputfile=str(out), extra_args=["--standalone"])
print(f"Created: {out}")
