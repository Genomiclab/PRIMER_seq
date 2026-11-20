# Input formats

## Cohort manifest

The batch drivers read a tab-separated manifest with one specimen per participant. Paths are absolute or relative to the manifest file.

| Column | Values |
| --- | --- |
| `sample_id` | Unique text key; leading zeros are preserved. Batch drivers accept letters, digits, underscores, dots and hyphens |
| `participant_id` | Participant key; include it when sample keys alone do not identify participants |
| `cohort` | `train`, `validation` or `blind` |
| `sex` | `F`, `M` or `U`; blank/NA is unknown |
| `group` | Training diagnosis: `Normal` or a supported cancer label |
| `mhap`, `cgi` | For `processed_cohorts.py`, fill exactly one per row |
| `bam`, `profile` | For `cohorts.py`, fill exactly one per row |
| `EM_valid_fraction`, `BPM_valid_fraction` | Required when canonical-only counts or canonical frequencies have discarded N-containing motifs |
| `methylation_block` | Optional training block-order key; fill every training row if the column is present |

Supported cancer labels, in model order: `Lung`, `Colorectum`, `Stomach`, `Liver`, `Thyroid`, `Breast`, `Lymphoma`, `Esophagus`, `Pancreas`, `Bladder`, `Prostate`, `Ovary`. The model includes only the classes present in training. Controls and at least two cancer types are required. Each group needs at least five samples for fold construction; meaningful training and tuning require more.

Validation and blind diagnoses belong in separate truth TSVs with `sample_id`, `group` and optional `stage` (`I`–`IV`). Prediction inputs do not need these labels. Step 5 requires exactly matching sample sets between predictions and truth.

## mHap and CGI statistics

mHap files must be nonempty, sorted, bgzip-compressed, six-column files aligned to hg19. Keep any existing tabix index adjacent. A missing index is created in the output workspace without modifying the source directory. The full hg19 CpG reference and its `.tbi` index are needed to summarize real samples.

To summarize one mHap file:

```bash
python3 scripts/step1_mhap_to_cgi.py \
  --mhap /path/to/sample.mhap.gz --cpg /path/to/hg19_CpG.gz \
  --out work/sample_cgi
```

The default BED is `config/hg19_cpgIsland.bed`. `cgi.tsv` contains `chr`, `start`, `end`, `nReads` and the methylation metrics. All samples must use the same region set. Step 2 retains autosomes and treats zero coverage as missing; model input replaces missing values with zero.

If CGI statistics already exist, use the `cgi` column in the manifest. No BAM, mHap conversion or genome package is then needed for this import.

## Cohort-wide EM and BPM matrices

Supply one matrix per channel. Text files may be TSV or CSV, optionally gzip-compressed. Their first column holds row keys and the header holds column keys. Matrices may have samples in rows or columns; the importer matches manifest sample keys to infer the orientation. Use `--layout samples-by-features` or `--layout features-by-samples` if needed. Extra samples are ignored; all requested samples must be present.

RDS files must contain one named numeric matrix or data frame. RData/RDA files must contain exactly one such object. A data frame may use a `sample_id` column instead of row names. Column and row names must be unique.

Include all 4,096 canonical A/C/G/T hexamers, including unobserved ones. Raw matrices can also contain N-containing hexamers; include all observed N-containing counts so their total is correct. Unobserved entries may be zero or NA. Negative values are rejected.

Choose the format explicitly when frequencies are supplied:

| `--motif-format` | Input values | Valid fraction |
| --- | --- | --- |
| `counts` (default) | Nonnegative integer counts | Computed from canonical and N-containing counts; required in metadata if only canonical counts remain |
| `canonical-frequencies` | Canonical counts divided by the canonical total; each row sums to 1 | Required in metadata |
| `all-frequencies` | Counts divided by the total including N-containing motifs | Computed from the canonical sum; N-containing columns may be present or omitted |

For each sample and channel, the valid fraction is `canonical_count / total_count_including_N`. Set it to 1 only when no ambiguous motifs were counted. Do not infer it from frequencies that have already been renormalized over canonical motifs.

Step 2 stores canonical frequencies and the valid fractions. Training and validation use the canonical frequencies; blind prediction multiplies by the valid fraction once, following the original analysis. All-frequency input is converted to this internal representation first.

To import one cohort with existing CGI tables, use a manifest with `sample_id`, `sex`, `cgi`, training `group` when applicable, and any required valid fractions:

```bash
Rscript scripts/step2_build_matrices.R \
  --cgi /path/to/cgi_manifest.tsv --em /path/to/EM.RData --bpm /path/to/BPM.RData \
  --motif_format counts --out work/matrices/train
```

The R command uses `--motif_format`; the Python batch command uses `--motif-format`.

## Seven precomputed feature tables

`step2_build_matrices.R --tables DIRECTORY` imports `samples.tsv` plus one `CHANNEL.tsv.gz` for each of MM, PDR, CHALM, MHL, Entropy, EM and BPM. Each table begins with `sample_id`; remaining columns are feature names. Row order is recovered by sample key.

For this input route, EM/BPM tables must contain all 4,096 canonical frequencies summing to one, and `samples.tsv` must include both valid fractions. Missing methylation values are allowed. CGI feature names must match the exact `chr.start.end` identifiers in the marker panel.

## Model inputs and outputs

Step 2 writes seven named numeric R matrices and `samples.tsv`. Step 3 reads `marker_file/markers.tsv` by default. Features absent from the matrices cause an error before fitting; measured regions with no coverage remain valid columns with missing values.

Step 3 writes new fold models, feature lists, fold assignments, out-of-fold probabilities, the training-derived threshold and `sessionInfo.txt`. Step 4 reads this output directory to score held-out samples. It rejects training/scoring sample overlap and, when supplied, participant overlap.

To inspect the imputed values supplied to a learner:

```bash
Rscript scripts/export_matrix.R work/matrices/blind blind EM work/blind_EM.tsv.gz
```
