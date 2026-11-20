# Installation and tested versions

Run on Linux or macOS with Python and R on `PATH`. The complete demonstration also requires the hg19 genome package for reference-based EM/BPM extraction. Use the same R installation for preprocessing, training and scoring. Each trained model records `sessionInfo.txt`.

## BAM tools

The tested command-line tools are mHapTools 0.10, mHapTk 1.0, samtools 1.3.1 and tabix 1.9. The wrapper uses Python's standard library. The historical mHapTk dependencies are pinned in `config/requirements-mhap.txt`.

```bash
conda create -n primer-mhap -c conda-forge -c bioconda python=3.7 pip samtools htslib
conda activate primer-mhap
python -m pip install -r config/requirements-mhap.txt
```

Install the C++ **mHapTools** program separately and put `mhaptools` on `PATH`. Its build needs a C++ compiler, HTSlib and compression libraries. The Python installation above supplies the separate `mhaptk` command. On Linux, building `pytabix` can require a C compiler and zlib development headers. [mHapTools build instructions](https://jiantaoshi.github.io/mHap/Installation.html) · [mHapTk source](https://github.com/CCCKW/mHapTk)

For full BAMs, obtain **both** `hg19_CpG.gz` and `hg19_CpG.gz.tbi` from the [upstream annotation page](https://jiantaoshi.github.io/mHap/AnnotationFiles.html). Keep them adjacent. The included `config/hg19_cpgIsland.bed` contains 27,718 intervals. The small reference in `examples/bam/` is sufficient only for its simulated BAM.

## R packages

The complete three-cancer demonstration was exercised with these versions:

| Component | Tested version |
| --- | --- |
| R | 4.1.3 |
| Bioconductor series | 3.14 |
| Rsamtools / GenomicAlignments | 2.10.0 / 1.30.0 |
| GenomeInfoDb / Biostrings | 1.30.1 / 2.62.0 |
| BSgenome.Hsapiens.UCSC.hg19 | 1.4.3 |
| data.table | 1.16.4 |
| caret / glmnet / gbm | 7.0-1 / 4.1-10 / 2.2.3 |
| e1071 | 1.7-17 |
| xgboost | 1.7.5.1 |
| catboost | 1.2 |
| doParallel / foreach | 1.0.17 / 1.5.2 |

Install in a dedicated R library. The following versions describe the tested environment; this is not a complete transitive dependency lock. Source installations require the compilers and headers for the packages concerned.

```r
# R 4.1.x
install.packages(c("BiocManager", "remotes"))
BiocManager::install(version = "3.14", ask = FALSE)
BiocManager::install(c("Rsamtools", "GenomicAlignments", "Biostrings",
  "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg19"), ask = FALSE)
versions <- c("data.table"="1.16.4", "e1071"="1.7-17", "gbm"="2.2.3",
  "glmnet"="4.1-10", "caret"="7.0-1", "xgboost"="1.7.5.1",
  "doParallel"="1.0.17", "foreach"="1.5.2")
for (p in names(versions))
  remotes::install_version(p, version=versions[[p]], upgrade="never")
# Install the CatBoost 1.2 R package for your operating system separately:
# install.packages("/path/to/catboost-R-package.tgz", repos=NULL)
```

CatBoost R packages are provided in the [official v1.2 release](https://github.com/catboost/catboost/releases/tag/v1.2). This implementation uses the XGBoost 1.7 R interface; its environment checker rejects incompatible major interfaces.

```bash
Rscript scripts/check_environment.R bam
Rscript scripts/check_environment.R train
Rscript scripts/check_environment.R all
```

The demonstration uses all five learners. Training time depends on the numerical libraries and CPU; the synthetic cohort is intentionally small. Full-size cohort matrices and models need substantially more memory and disk space. Use `--threads` for learner parallelism and `--jobs` for independent BAM processing.
