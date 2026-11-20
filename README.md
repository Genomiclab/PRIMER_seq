# PRIMER_seq: BAM to ensemble training and evaluation

This repository provides a step-by-step workflow for paired-end PRIMER plasma DNA data: BAM → mHap, end motifs (EM), breakpoint motifs (BPM) → five methylation channels and two fragmentomic channels → two-layer ensemble training → validation and blinded prediction.

The release includes a **simulated BAM**, **synthetic Lung, Colorectum, Liver and Normal cohorts**, and a small marker panel. All model weights are fitted when the example runs. **The original study's fitted models and participant-level data are not distributed in this repository.** The examples demonstrate the computation; their results do not reproduce or estimate the manuscript's clinical performance.

[中文说明](docs/QUICKSTART_zh.md) · [Installation](docs/ENVIRONMENT.md) · [Inputs and methods](docs/METHODS.md) · [Example provenance and checks](docs/VERIFICATION.md)

## Run the complete example

Install the dependencies, then run from this repository's root:

```bash
Rscript scripts/check_environment.R all
python3 examples/run_demo.py --out work/demo --threads 1
```

The driver runs all five steps, trains all five algorithms in all five folds, predicts the BAM-derived sample, and evaluates the separate validation and blind cohorts. Every command requires a new output directory. Outputs are kept under `work/`, which is excluded from Git.

| Example | Contents |
| --- | --- |
| `examples/bam/` | 3,880 simulated read pairs across 97 hg19 CpG islands, BAM index, small CpG annotation and BED |
| `examples/three_cancers/train/` | 160 synthetic samples: 40 per group |
| `examples/three_cancers/validation/` | 40 synthetic samples: 10 per group |
| `examples/three_cancers/blind/` | 40 synthetic samples: 10 per group |
| `examples/three_cancers/markers.tsv` | Supplied feature identities for detection and tissue of origin |
| `examples/expected_step1/` | Reference outputs for the supplied simulated BAM |
| `examples/expected_results/` | Aggregate example results from a tested environment |

The single BAM illustrates sample processing. The separate synthetic cohorts provide enough independent rows for training and evaluation. For an actual cohort, process each participant's BAM using Step 1 and assemble those profiles using Step 2.

## Step 1 — Generate mHap, methylation features, EM and BPM

```bash
python3 scripts/step1_bam_to_features.py \
  --bam examples/bam/synthetic.bam \
  --cpg examples/bam/hg19_demo_CpG.gz \
  --bed examples/bam/regions.bed \
  --out work/step1
```

Outputs are indexed `sample.mhap.gz`, `cgi.tsv`, combined raw counts in `motifs.tsv`, and separate `EM.tsv` and `BPM.tsv` files. The five methylation channels used by the model are **MM, PDR, CHALM, MHL and Entropy**. Each fragmentomic channel uses 4,096 A/C/G/T hexamers.

For your samples, supply coordinate-sorted, paired-end, methylation-informative PRIMER BAMs aligned to hg19, the full indexed hg19 CpG reference, and `config/hg19_cpgIsland.bed`. The small example reference covers only the demonstration regions. [Reference and software setup](docs/ENVIRONMENT.md)

## Step 2 — Assemble sample-by-feature matrices

Assemble the profile generated in Step 1:

```bash
Rscript scripts/step2_build_matrices.R \
  --profiles examples/bam_sample.tsv --out work/matrices/bam
```

Import the supplied three-cancer training, validation and blind tables:

```bash
for cohort in train validation blind; do
  Rscript scripts/step2_build_matrices.R \
    --tables examples/three_cancers/$cohort --out work/matrices/$cohort
done
```

Each cohort directory contains seven numeric matrices and `samples.tsv`. Training metadata include `group`. Validation and blind diagnoses are stored separately in `truth.tsv` and are read only by Step 5. [Metadata and matrix schema](docs/METHODS.md)

## Step 3 — Train the two-layer ensemble with supplied markers

```bash
Rscript scripts/step3_train_ensemble.R \
  --train work/matrices/train \
  --markers examples/three_cancers/markers.tsv \
  --out work/model --preset demo --threads 1
```

The first layer trains **elastic-net logistic regression, GBM, linear SVM, XGBoost and CatBoost**, separately for each channel and task. The second layer trains a linear SVM on first-layer out-of-fold probabilities. Detection and tissue-of-origin models each use five folds; the complete example produces 360 newly fitted fold models.

The `demo` preset reduces internal tuning and boosting iterations for this small dataset. The `study` preset uses the original script settings, including 1,000 XGBoost rounds and 2,000 CatBoost iterations. Both presets use all five algorithms and the same two-layer architecture. To train on your cohort, supply its matrices and marker table and use `--preset study`. [Exact settings](docs/METHODS.md)

The output contains newly trained weights, selected feature lists, fold assignments, out-of-fold predictions, a newly estimated detection threshold, and the R session information.

## Step 4 — Predict validation and blind samples

```bash
Rscript scripts/step4_predict.R \
  --data work/matrices/validation --model work/model \
  --cohort validation --out work/predictions/validation

Rscript scripts/step4_predict.R \
  --data work/matrices/blind --model work/model \
  --cohort blind --out work/predictions/blind

Rscript scripts/step4_predict.R \
  --data work/matrices/bam --model work/model \
  --cohort validation --out work/predictions/bam
```

Predictions contain detection scores and calls, the fixed threshold learned in Step 3, and tissue probabilities for detection-positive samples. The third command completes the BAM-to-prediction example; that simulated BAM has no clinical diagnosis. Inference never reads `truth.tsv`.

## Step 5 — Evaluate after predictions are fixed

```bash
for cohort in validation blind; do
  Rscript scripts/step5_evaluate.R \
    --predictions work/predictions/$cohort/predictions.tsv \
    --truth examples/three_cancers/$cohort/truth.tsv \
    --out work/reports/$cohort
done
```

Reports include confusion counts, sensitivity, specificity, AUC, ROC coordinates, tissue-of-origin top-k accuracy and confusion counts, and a PDF plot. Tissue accuracy is calculated among detected true cancer cases.

## Process your own cohort

Fill [the manifest template](examples/cohorts_template.tsv) with one specimen per participant and existing train/validation/blind assignments. The template rows are illustrative; they are not a training dataset. Supply either a BAM or an already completed Step 1 profile per row.

```bash
python3 workflow/cohorts.py --manifest /path/to/cohorts.tsv \
  --cpg /path/to/hg19_CpG.gz --bed config/hg19_cpgIsland.bed \
  --out work/cohorts --jobs 2

Rscript scripts/step3_train_ensemble.R \
  --train work/cohorts/matrices/train --markers /path/to/markers.tsv \
  --out work/cohort_model --preset study --threads 4
```

Continue with Steps 4–5 using `work/cohorts/matrices/validation` and `work/cohorts/matrices/blind`. Truth labels remain outside the inference inputs. The batch entry checks sample and participant uniqueness across all three cohorts. Large cohorts can run Step 1 as separate cluster jobs and reuse their completed profiles.

Run `Rscript tests/motif_denominators.R` and `Rscript tests/participant_separation.R` for focused input and scoring checks. To regenerate the simulated BAM, run `Rscript examples/make_demo_bam.R work/new_simulated_bam`.
