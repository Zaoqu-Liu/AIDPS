# AIDPS

Artificial intelligence-derived prognostic signature for pancreatic cancer (PACA).

## Citation

Wang L, Liu Z, Liang R, Wang W, Zhu R, Li J, Xing Z, Weng S, Han X, Sun Y. Comprehensive machine-learning survival framework develops a consensus model in large-scale multicenter cohorts for pancreatic cancer. *eLife* 2022;11:e80150. doi: [10.7554/eLife.80150](https://doi.org/10.7554/eLife.80150)

## Summary

As the most aggressive tumor, the outcome of pancreatic cancer (PACA) has not improved observably over the last decade. Anatomy- based TNM staging does not exactly identify treatment- sensitive patients, and an ideal biomarker is urgently needed for precision medicine. Based on expression files of 1280 patients from 10 multicenter cohorts, we screened 32 consensus prognostic genes. Ten machine- learning algorithms were transformed into 76 combinations, of which we selected the optimal algorithm to construct an artificial intelligence- derived prognostic signature (AIDPS) according to the average C- index in the nine testing cohorts. The results of the training cohort, nine testing cohorts, Meta- Cohort, and three external validation cohorts (290 patients) consistently indicated that AIDPS could accurately predict the prognosis of PACA. After incorporating several vital clinicopathological features and 86 published signatures, AIDPS exhibited robust and dramatically superior predictive capability. Moreover, in other prevalent diges-tive system tumors, the nine- gene AIDPS could still accurately stratify the prognosis. Of note, our AIDPS had important clinical implications for PACA, and patients with low AIDPS owned a dismal prognosis, higher genomic alterations, and denser immune cell infiltrates as well as were more sensitive to immunotherapy. Meanwhile, the high AIDPS group possessed observably prolonged survival, and panobinostat may be a potential agent for patients with high AIDPS. Overall, our study provides an attractive tool to further guide the clinical management and individualized treat-ment of PACA.

## Repository contents

| File | Purpose |
| --- | --- |
| `AIDPS.R` | Builds the signature: CoxBoost selects genes, then survival SVM (`survivalsvm`) fits the risk score on the training cohort and predicts the validation cohort |
| `AIDPS.rda` | Result of `AIDPS.R`: `rid` (the 9 selected genes, Ensembl IDs) and `fit` (the fitted `survivalsvm` model) |
| `PACA_AU_Array.rda` | Training cohort, a data frame of 267 patients: `ID`, `OS.time`, `OS`, and 32 genes (Ensembl IDs) |
| `TCGA_PAAD.rda` | Validation cohort, 176 patients, same columns |

## Usage

Requires R with `survival` and `survivalsvm`; re-running `AIDPS.R` also needs `CoxBoost` and `tidyverse`. Run from the repository directory.

Score patients with the saved model:

```r
library(survival)
library(survivalsvm)

load('AIDPS.rda')      # rid (selected genes), fit (survivalsvm model)
load('TCGA_PAAD.rda')  # or your own data frame with the genes in `rid` as columns

RS <- as.numeric(predict(fit, TCGA_PAAD[, rid])$predicted)
```

A higher score corresponds to longer survival (on `TCGA_PAAD`, Cox regression of overall survival on the score gives a coefficient below zero, p = 6e-06), which is why a low AIDPS marks the poor-prognosis group. Columns must be expression values of the genes in `rid`, named by the same Ensembl IDs and on the same scale as the training data.

Rebuild the model from the training cohort:

```r
source('AIDPS.R')  # reads PACA_AU_Array.rda and TCGA_PAAD.rda, writes AIDPS.rda
```
