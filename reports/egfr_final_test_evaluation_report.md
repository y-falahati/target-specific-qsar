Final Test Evaluation --- EGFR Bioactivity Prediction
Project
Project: Target-Specific Bioactivity Prediction and Virtual
Screening Using ChEMBL, RDKit and Machine Learning
Target: Human EGFR (CHEMBL203; UniProt P00533)
Evaluation stage: Final held-out Test evaluation
Test set size: 697 compounds
Objective
This evaluation estimates how well the saved Random Forest models
predict compound pIC50 values for compounds held out from model
training, hyperparameter tuning, and validation-based model selection.
Evaluation protocol
- The pre-defined scaffold-based split was used.
- The Test subset was used for prediction and final reporting only.
- The saved baseline and tuned Random Forest models were loaded
  without retraining.
- No hyperparameters were changed during this evaluation.
- Metrics were calculated from observed and predicted pIC50 values for
  the 697 Test compounds.
The tuned model had previously been selected as a candidate based on the
training/validation workflow. The baseline model is included here as a
reference comparison.
Final Test results
  Metric                                   Random Forest baseline   Random Forest tuned
  Number of Test compounds                                    697                   697
  R²                                                     0.250845              0.284709
  RMSE                                                   1.131138              1.105277
  MAE                                                    0.862515              0.859525
  Mean residual (observed − predicted)                  -0.184470             -0.190407
  Median absolute error                                  0.664938              0.702812
  Maximum absolute error                                 3.960670              3.955282
Interpretation
The tuned Random Forest achieved an R² of 0.2847, RMSE of
1.1053, and MAE of 0.8595 on the held-out Test set.

Compared with the baseline model:
- R² increased from 0.250845 to 0.284709.
- RMSE decreased from 1.131138 to 1.105277, an approximately 2.29%
  reduction.
- MAE decreased from 0.862515 to 0.859525, an approximately 0.35%
  reduction.
- Median absolute error increased from 0.664938 to 0.702812, so the
  tuned model did not improve every error summary.
- The mean residual was negative for both models. With residual
  defined as observed pIC50 minus predicted pIC50, this indicates
  average overprediction of pIC50.
Overall, tuning produced a modest improvement in the principal
aggregate metrics, rather than a large performance gain. The R²
indicates limited explanatory predictive performance on this Test set.
The results should not be described as high-accuracy bioactivity
prediction.

Because pIC50 is logarithmic, an error of one pIC50 unit corresponds to
a tenfold difference in the associated concentration scale. The observed
error magnitudes are therefore meaningful for medicinal-chemistry
interpretation.
Limitations
1. This is a held-out Test evaluation for the current split and the two
   saved models; it is not an external validation on an independent
   dataset.
2. The modest difference between models has not been established as
   statistically significant.
3. The model uses 13 RDKit molecular descriptors, which provide a
   compact representation and may not capture all structural
   determinants of EGFR activity.
4. The curated activity labels aggregate measurements across available
   records at the compound level; experimental heterogeneity may
   contribute to prediction error.
5. The negative mean residual indicates average overprediction, but
   does not establish that the model is biased in the same direction
   across every activity range.
Reproducibility and output files
Evaluation script:

notebooks/18_final_test_evaluation.py

Metrics:

results/tables/egfr_final_test_metrics.csv

Compound-level predictions:

results/tables/egfr_final_test_predictions.csv
Final conclusion
On 697 held-out Test compounds, the tuned Random Forest achieved R² =
0.2847, RMSE = 1.1053, and MAE = 0.8595. It modestly outperformed the
baseline on R², RMSE, and MAE, but had a higher median absolute error.
The model is a useful baseline for this target-specific QSAR workflow,
but its predictive accuracy remains limited and should be improved and
independently validated before use in decision-critical virtual
screening.

Test-set integrity note: These results are the final reported
results for this split. Do not use this Test set to tune features,
hyperparameters, or model choice and then continue to describe the same
set as an untouched final test set