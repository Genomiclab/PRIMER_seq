# PRIMER_seq

Processing and ensemble modelling of paired-end PRIMER plasma DNA data. The workflow uses five methylation features (MM, PDR, CHALM, MHL and Entropy), end motifs (EM) and breakpoint motifs (BPM) for cancer detection and tissue-of-origin prediction.

The study marker lists are in [`marker_file/`](marker_file/README.md). Training uses these lists by default. Example samples are simulated; participant data and the original fitted models are not included.

[Installation](docs/ENVIRONMENT.md) · [Input formats](docs/INPUTS.md) · [Methods](docs/METHODS.md) · [中文说明](docs/QUICKSTART_zh.md)

## Run the example

After installing the dependencies, run from the repository root:

```bash
Rscript scripts/check_environment.R all
python3 examples/run_demo.py --out work/demo --threads 1
```

This runs Steps 1–5 with the complete marker panel, all five learners and five folds. It processes one simulated BAM, trains on 160 synthetic samples (Normal, Lung, Colorectum and Liver), and evaluates separate validation and blind sets of 40 samples each. The BAM-derived sample is also scored. All model weights are fitted during the run.

The supplied small synthetic tables are expanded to the full methylation marker panel by `examples/make_study_demo.R`. Generated values are simulated, including values at the study marker coordinates. Full-panel training is CPU intensive: 360 fold models are fitted even with the shortened demo preset. Use a new output directory for each run.

## Step 1: BAM to mHap, methylation features, EM and BPM

```bash
python3 scripts/step1_bam_to_features.py \
  --bam examples/bam/synthetic.bam \
  --cpg examples/bam/hg19_demo_CpG.gz \
  --bed config/hg19_cpgIsland.bed \
  --out work/step1
```

The output directory contains `sample.mhap.gz` and its index, `cgi.tsv`, raw counts in `motifs.tsv`, and separate `EM.tsv` and `BPM.tsv` files. The model uses autosomal CGI features and canonical six-base motifs.

For other BAMs, supply the full indexed hg19 CpG reference. The small reference shipped with the example covers only the simulated reads. Use the supplied full CGI BED to retain every marker region.

## Step 2: Build sample-by-feature matrices

Assemble the single BAM profile:

```bash
Rscript scripts/step2_build_matrices.R \
  --profiles examples/bam_sample.tsv --out work/matrices/bam
```

Prepare and import the synthetic cohorts:

```bash
Rscript examples/make_study_demo.R --out work/tables
for cohort in train validation blind; do
  Rscript scripts/step2_build_matrices.R \
    --tables work/tables/$cohort --out work/matrices/$cohort
done
```

Each output contains seven R matrices and `samples.tsv`. Training metadata include `group`; held-out diagnoses are kept in separate `truth.tsv` files.

## Step 3: Train the ensemble

```bash
Rscript scripts/step3_train_ensemble.R \
  --train work/matrices/train --out work/model --preset demo --threads 1
```

The default panel is `marker_file/markers.tsv`, with separate detection and tissue-of-origin lists for all seven channels. The script checks that every marker is present before fitting. The ordered features actually used are saved as `markers_used.tsv` in the model directory.

The first layer fits elastic-net logistic regression, GBM, linear SVM, XGBoost and CatBoost. A second-layer linear SVM combines their out-of-fold probabilities. Each task uses five folds at each layer. The `demo` preset reduces tuning and boosting iterations; `--preset study` uses the original script settings. See [training settings](docs/METHODS.md#ensemble-and-training-presets).

## Step 4: Predict validation and blind samples

```bash
for cohort in validation blind; do
  Rscript scripts/step4_predict.R \
    --data work/matrices/$cohort --model work/model \
    --cohort $cohort --out work/predictions/$cohort
done
```

Predictions include detection scores, calls at the training-derived threshold, and tissue probabilities for detection-positive samples. `--cohort blind` applies the original blind-cohort motif normalization. Inference does not read truth labels.

The single BAM can be scored with the same command using `--data work/matrices/bam --cohort validation --out work/predictions/bam`.

## Step 5: Evaluate fixed predictions

```bash
for cohort in validation blind; do
  Rscript scripts/step5_evaluate.R \
    --predictions work/predictions/$cohort/predictions.tsv \
    --truth work/tables/$cohort/truth.tsv --out work/reports/$cohort
done
```

Reports contain sensitivity, specificity, AUC, ROC coordinates, tissue top-k accuracy, confusion counts and a PDF. Synthetic example results are software checks, not clinical performance estimates.

## Start from existing mHap and EM/BPM matrices

Fill [`examples/processed_cohorts_template.tsv`](examples/processed_cohorts_template.tsv) with one specimen per participant and the existing cohort assignments. Use one mHap file or CGI statistics file per row. EM/BPM matrices may be TSV, CSV, RDS or single-object RData files; see [formats and normalization](docs/INPUTS.md).

```bash
python3 workflow/processed_cohorts.py \
  --manifest /path/to/cohorts.tsv \
  --em /path/to/EM.tsv.gz --bpm /path/to/BPM.tsv.gz \
  --motif-format counts --cpg /path/to/hg19_CpG.gz \
  --out work/cohorts --jobs 2

Rscript scripts/step3_train_ensemble.R \
  --train work/cohorts/matrices/train --out work/cohort_model \
  --preset study --threads 4
```

Continue with Steps 4–5 using `work/cohorts/matrices/validation` and `work/cohorts/matrices/blind`, the newly fitted `work/cohort_model`, and the corresponding truth tables. The batch driver checks sample and participant separation across cohorts. If starting from BAMs, use [`workflow/cohorts.py`](workflow/cohorts.py) with [`examples/cohorts_template.tsv`](examples/cohorts_template.tsv).

Input checks and reference outputs are described in [verification](docs/VERIFICATION.md). Generated matrices, models and predictions belong under `work/`.
