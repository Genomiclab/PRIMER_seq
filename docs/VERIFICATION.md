# Synthetic examples and checks

## Example inputs

`examples/bam/` contains 3,880 simulated read pairs across 97 hg19 CpG islands. The BAM is generated from the public reference using simulated fragment positions, methylation states and bisulfite-converted bases. Read names begin `SIMFRAG`. Regenerate it with:

```bash
Rscript examples/make_demo_bam.R work/new_simulated_bam
```

`examples/three_cancers/` contains fixed synthetic seed tables for Normal, Lung, Colorectum and Liver. Training has 160 samples, validation 40 and blind testing 40, with equal group sizes. Identifiers, sex and stage labels are simulated. Held-out diagnoses are stored in separate truth files.

`examples/make_study_demo.R` expands the five methylation tables to the complete study marker panel. It uses seeded mixtures of existing synthetic columns and added noise, bounded to [0,1]. The same column mixtures are used across cohorts; noise is generated separately. EM/BPM already contain all 4,096 canonical motifs. This generator reads only public example tables and marker identities.

The small seed panel in `examples/three_cancers/markers.tsv` is used to construct the initial example, including the 97-region BAM. Default model training uses `marker_file/markers.tsv`.

## Reference outputs

`examples/expected_step1/` contains the simulated BAM outputs with the full `config/hg19_cpgIsland.bed`. Most regions have no coverage in this small BAM. Compare uncompressed mHap records and numerical contents rather than gzip headers.

The tested software versions are in [ENVIRONMENT.md](ENVIRONMENT.md). Full-panel checks covered marker coverage, matrix assembly, binary and multiclass learner fitting, and SVM interface equivalence. The complete 360-model training and held-out scoring run has not been completed for this marker update. Steps 4–5 generate predictions and evaluation tables when the demonstration is run.

## Input checks

```bash
Rscript tests/motif_denominators.R
Rscript tests/participant_separation.R
Rscript tests/processed_inputs.R
Rscript tests/svm_input.R
```

These checks cover canonical/all-count denominators, missing motifs and zero coverage, text sample keys, sex handling, participant overlap, matrix orientation, raw/frequency imports and equivalence with the BAM-profile input route.

The mHap batch route was also checked with synthetic inputs, including an unindexed source mHap file, existing CGI statistics, reordered motif matrices and separate cohort manifests. The source file was preserved, and held-out diagnoses were omitted from the resulting prediction metadata.

To check the merged study panel against the source CSVs:

```bash
python3 scripts/build_markers.py --out work/marker_check
diff marker_file/markers.tsv work/marker_check/markers.tsv
diff marker_file/counts.tsv work/marker_check/counts.tsv
```

Training checks the full selected panel before fitting. A missing marker column stops the run without writing a model directory. Step 3 records the ordered `markers_used.tsv`, fold assignments and R session information on completion.
