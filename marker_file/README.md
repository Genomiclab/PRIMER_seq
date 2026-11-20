# Study markers

These 14 CSV files contain the feature identities used for cancer detection and tissue of origin. The two tasks have separate lists. No participant measurements or fitted model parameters are stored here.

| Channel | Cancer detection | Tissue of origin |
| --- | ---: | ---: |
| MM | 11,841 | 10,381 |
| PDR | 12,953 | 11,331 |
| CHALM | 10,133 | 9,094 |
| MHL | 11,028 | 9,566 |
| Entropy | 13,233 | 10,686 |
| EM | 3,999 | 3,983 |
| BPM | 3,978 | 3,992 |

Methylation CSVs contain `metric` and `hg19_CGI_genomic_region`; regions use `chr.start.end` identifiers matching `config/hg19_cpgIsland.bed`. Motif CSVs contain `metric` and `motif`; motifs are six-base A/C/G/T sequences. Any leading CSV row-number column is ignored.

`markers.tsv` combines the CSVs into `task`, `channel`, `feature` columns and preserves the feature order within each source file. `task` is `detection` or `origin`. This is the default input to Step 3, including the synthetic demonstration. `counts.tsv` records the source filename and count for each list.

To rebuild the combined tables in a new directory:

```bash
python3 scripts/build_markers.py --input marker_file --out work/rebuilt_markers
```

Step 3 accepts `--markers /path/to/markers.tsv` for a different panel. Every task/channel list must have at least two distinct features, and every selected feature must be present in the input matrices.
