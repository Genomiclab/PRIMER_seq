# Example provenance and verification

## Public example data

The BAM is generated from the public hg19 reference using simulated fragment positions, methylation states and bisulfite-converted bases. It contains 3,880 read pairs, 85-base reads, and 97 demonstration CpG islands. Read names begin `SIMFRAG`; it contains no sampled participant reads or original sample identifiers. `examples/make_demo_bam.R` and `examples/simulation.tsv` regenerate this example with seed 20261120. Read/fragment-length settings were chosen from rounded summaries of a private plasma check.

The three-cancer numerical tables are independently drawn synthetic observations. Their ranges were calibrated from feature means and variances in the original training subset for Normal, Lung, Colorectum and Liver. Noise and shared simulated latent terms were added, methylation values were bounded to [0,1], and canonical motif frequencies were normalized to sum to one. No original participant row was copied. Sex and stage values are artificial. The aggregate calibration inputs remain private; the supplied synthetic tables are the fixed public inputs for the example.

Training, validation and blind examples contain distinct synthetic sample keys. Held-out diagnosis and stage labels reside in `truth.tsv`; the prediction command receives no truth file. Because the example labels are openly available, this demonstrates the blind-scoring interface rather than real clinical blinding.

## Completed checks

The public `examples/run_demo.py` command completed from the supplied inputs, including:

- BAM → indexed mHap, methylation region statistics, EM and BPM.
- Step 2 assembly of the BAM profile and import of all three synthetic cohorts.
- All five algorithms, all five folds, both tasks and both model layers.
- Prediction of the BAM-derived synthetic sample using the newly fitted model.
- Separate validation and blind predictions, followed by truth-based evaluation.

The example output has three trained tissue classes and 360 newly fitted fold files. Fixed-input results in the tested single-thread environment were:

| Synthetic cohort | n | TP | FN | TN | FP | AUC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Validation | 40 | 20 | 10 | 10 | 0 | 0.9300 |
| Blind | 40 | 25 | 5 | 10 | 0 | 0.9733 |

These are software demonstration results, not study performance estimates. The complete aggregate outputs are under `examples/expected_results/`; exact probabilities can vary across environments. `examples/expected_step1/` contains reference outputs for the simulated BAM. Compare numerical contents, since compressed file headers can vary.

Focused regression checks cover N-containing motif denominators, zero coverage, missing motifs, text sample keys, sex adjustment with reduced class sets, and rejection of training/scoring participant overlap.

The source BAM processing was also exercised on a private 1,500-fragment plasma subset on macOS and a Linux HPC environment. The uncompressed mHap records and all EM/BPM count and frequency tables matched exactly. The five modelling metrics matched for all covered regions in a 40-region check; zero-coverage regions are consistently converted to missing by Step 2. The installed mHapTk versions differed in auxiliary output columns and zero-coverage text representation. That subset and its derived files are excluded from this release. The original study scripts remain separate from the public implementation.

## Release contents

The repository contains source code, configuration, a reduced marker panel, simulated inputs and their example outputs. It contains no original fitted model, serialized R model object, original control threshold, original sample list or clinical table. Runtime outputs belong under the Git-ignored `work/` directory. Model output directories contain fitted parameters and should be handled according to the user's own data-release policy.
