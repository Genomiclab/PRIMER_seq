# Inputs and computational methods

## Sample and feature tables

| Field | Meaning |
| --- | --- |
| `sample_id` | Unique sample key, preserved as text, including leading zeros |
| `participant_id` | Optional participant key; provide it when different sample keys could belong to the same participant |
| `sex` | `F`, `M` or `U` for unknown; blank/NA is treated as unknown by the R reader |
| `group` | Required training diagnosis; evaluation diagnoses belong in the separate truth table |
| `stage` | Optional truth-table field: `I`, `II`, `III`, `IV`, or missing |
| `profile` | Directory from Step 1, absolute or relative to the manifest TSV |
| `mhap`, `cgi` | Alternative per-sample inputs for `workflow/processed_cohorts.py`; see [input formats](INPUTS.md) |
| `cohort` | Batch-driver assignment: `train`, `validation`, `blind` |
| `bam` | Batch-driver input BAM; supply exactly one of `bam` and `profile` per row |
| `EM_valid_fraction`, `BPM_valid_fraction` | Fraction of total raw motif counts with no N, retained for the original blind-cohort normalization |
| `methylation_block` | Optional original training batch-order key, used for ordering methylation rows before fold generation |

Supported cancer labels, in model order, are `Lung`, `Colorectum`, `Stomach`, `Liver`, `Thyroid`, `Breast`, `Lymphoma`, `Esophagus`, `Pancreas`, `Bladder`, `Prostate`, `Ovary`. `Normal` is the non-cancer control label. A fit includes only the cancer classes present in its training cohort. At least two cancer types and controls are required. Each group needs at least five samples for the outer five-fold split; substantially more are needed for meaningful fitting and internal tuning. The example provides 40 per group.

Use one specimen per participant. The batch driver checks participant uniqueness across all cohorts and forwards participant keys to the matrices. Scoring rejects overlap with training sample keys, and with training participant keys when supplied. Keep validation and blind participants distinct as well.

`--profiles` reads each sample's `cgi.tsv` and `motifs.tsv`. `--tables` imports `samples.tsv` plus `MM.tsv.gz`, `PDR.tsv.gz`, `CHALM.tsv.gz`, `MHL.tsv.gz`, `Entropy.tsv.gz`, `EM.tsv.gz`, and `BPM.tsv.gz`. Each numerical table starts with a `sample_id` column, followed by named feature columns. Row order is recovered by sample key. EM/BPM tables must contain all 4,096 distinct A/C/G/T hexamers, with nonnegative frequencies summing to one per sample; supply both valid-fraction fields in the metadata. Missing methylation values are allowed.

Step 2 writes seven numeric R matrices with sample and feature names and a metadata TSV. To export the actual imputed input values for a specified cohort:

```bash
Rscript scripts/export_matrix.R work/matrices/blind blind EM work/blind_EM.tsv.gz
```

## Methylation and genomic coordinates

Inputs are coordinate-sorted, paired-end, methylation-informative PRIMER alignments against UCSC hg19. Preserve the original sequence and mate/strand flags. mHapTools `convert` reads the BAM and indexed CpG annotation, and mHapTk `stat` summarizes the requested BED regions. mHap records contain chromosome, first and last covered CpG position, haplotype string, count and strand. [mHapTools documentation](https://jiantaoshi.github.io/mHap/Installation.html) · [mHapTk source and commands](https://github.com/CCCKW/mHapTk)

The region table contains MM, PDR, CHALM, MHL, MCR, MBS and Entropy. The ensemble uses **MM, PDR, CHALM, MHL and Entropy** from chr1–chr22. Region identifiers are the exact `chr.start.end` coordinates from the BED/stat output. Zero-coverage regions remain missing until model input; missing values are then imputed to zero, following the source workflow. All samples in a profile assembly must have the same region set.

Step 1 uses a temporary working directory and creates a local index when the input BAM has no adjacent index. It does not change the input BAM or write into the original BAM directory.

## EM and BPM

The original calculation counts eligible mapped, nonsecondary, properly paired R2 alignments on chr1–chr22. It uses the hg19 reference span implied by the CIGAR. The source flags add no separate duplicate or supplementary-alignment filter. The MspI exclusion checks the R2 end; it does not require an additional positive CCGG test on R1. Coordinates below are one-based and inclusive. Minus-strand 6-mers are reverse-complemented.

| R2 position | EM interval | BPM interval | Reject if this interval is CCGG |
| --- | --- | --- | --- |
| Plus-strand start `s` | `s−6 … s−1` | `s−9 … s−4` | `s−7 … s−4` |
| Minus-strand end `e` | `e+1 … e+6` | `e+4 … e+9` | `e+4 … e+7` |

`motifs.tsv` contains all 15,625 A/C/G/T/N hexamers and raw counts for both channels. `EM.tsv` and `BPM.tsv` provide `kmer`, `count`, and canonical `frequency` separately. N-containing motifs have `NA` frequency. Canonical frequencies divide by the sum of counts without N.

The two model matrices retain the 4,096 canonical hexamers. Absent motifs are represented as missing, then imputed to zero at model input. Training and internal validation normalize by the canonical total. The original blinded script normalizes by all counts, including N-containing motifs. Step 2 records the valid fractions, and `--cohort blind` multiplies the canonical frequencies by them once. The supplied synthetic cohort has valid fractions of one; `tests/motif_denominators.R` exercises non-unit fractions.

## Study markers

The required marker TSV contains `task`, `channel`, `feature`. Tasks are `detection` and `origin`; channels are the seven names above. Preserve the intended order within a task/channel. Every list must contain at least two distinct features available in the corresponding matrix. No fitted coefficients or tree parameters are present in this file.

`marker_file/markers.tsv` is the default panel for training and the demonstration. It combines the 14 study CSVs, preserving separate detection and origin lists and their feature order. Counts and file formats are listed in [marker_file/README.md](../marker_file/README.md). Step 3 checks coverage of the complete selected panel before fitting and records the features in `markers_used.tsv`.

The small panel under `examples/three_cancers/markers.tsv` defines the seed synthetic tables and the 97-region BAM example. It is not the default training panel. `examples/make_study_demo.R` expands the public synthetic methylation tables to every region required by the study panel, using deterministic mixtures of seed columns and added noise. EM/BPM tables already cover all 4,096 canonical motifs. The generated values remain synthetic.

By default all supplied marker identities are fixed before fitting. For optional training-only motif selection, use `--select_motifs true`. In that mode, EM/BPM columns with fewer than 50 missing rows enter Wilcoxon comparisons, with BH-adjusted P < 0.01 within each comparison. Detection compares each cancer to controls; origin compares cancer pairs. Selection takes place within the complete training cohort before its folds are constructed, as in the source scripts. Validation and blind labels never enter feature selection or training. Use an externally specified marker panel or nested selection if a different validation design requires it.

## Ensemble and training presets

Each task has seven input channels, five first-layer learners and five folds: 175 first-layer fits per task. The second layer uses a separate five-fold linear SVM for each task. The complete run produces 350 first-layer fits and 10 second-layer fits. The tissue-of-origin task uses cancer cases only. First-layer out-of-fold probabilities form 35 detection inputs and `35 × number_of_cancer_types` origin inputs (105 for the example).

| Learner | `study` preset | `demo` preset |
| --- | --- | --- |
| LR | glmnet, alpha 0.5, 10 internal folds, `lambda.min` | Same, 5 internal folds |
| GBM | caret `gbm`, default tuning/resampling | Fixed 25 trees, depth 1, shrinkage 0.1, minimum node size 5; no inner resampling |
| Linear SVM | e1071, scaled features, probability estimates | Same |
| XGBoost | 1,000 rounds | 20 rounds |
| CatBoost | 2,000 iterations | 25 iterations |

First-layer seeds are 35 for detection and 21 for origin, reset per learner. The second-layer seeds are 2 and 100. Held-out predictions average the five fitted fold probabilities at each layer. R version, package versions, numerical libraries and thread counts can affect numerical results.

The linear SVM uses e1071's numeric matrix interface to avoid constructing a large formula term matrix for the full marker panel. With the same numeric inputs, labels and seed, the two interfaces produced identical support vectors, coefficients and probabilities in the binary and three-class checks.

The detection threshold is estimated anew from training-control out-of-fold scores at the 97th percentile using `quantile(type=1)`. Validation uses `score > threshold`; blinded scoring uses `score >= threshold`, preserving the source script behavior. This percentile choice does not guarantee 97% specificity in a separate cohort.

For origin, the source scoring rule zeros Ovary/Breast for samples coded M and Prostate for samples coded F, when those classes are present, and renormalizes. Unknown sex leaves the class probabilities intact. The three-cancer demonstration has none of those sex-specific columns. Origin predictions are generated only for detection-positive held-out samples.

## Evaluation and reproducibility scope

Step 5 joins prediction and truth files by unique sample key and requires exactly matching sample sets. Detection reports confusion counts, sensitivity, specificity, rank-based AUC and ROC coordinates, with optional stage summaries. Origin top-k accuracy is conditioned on detected true cancers, with k capped at the number of trained tissue classes.

The public inputs and training code reproduce the synthetic demonstration using the study marker identities. Participant measurements, cohort assignments and original fitted parameters are not included. A new fit produces new model weights and a training-derived threshold; numerical results depend on the input cohort and software environment.
