#!/usr/bin/env python3
"""Extract study features from one coordinate-sorted, paired-end hg19 PRIMER BAM."""
import argparse
import csv
import gzip
from itertools import zip_longest
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def preprocess(bam, cpg, bed, out):
    bam, cpg, bed = [p.resolve(strict=True) for p in (bam, cpg, bed)]
    out = out.absolute()
    if out.exists():
        raise ValueError(f"Use a new output directory: {out}")
    if not Path(str(cpg) + ".tbi").is_file():
        raise ValueError("Supply the bgzip-compressed CpG reference together with its .tbi index")
    for tool in ("mhaptools", "mhaptk", "tabix", "samtools", "Rscript"):
        if not shutil.which(tool):
            raise ValueError(f"{tool} must be on PATH")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".primer-", dir=out.parent) as tmp:
        work = Path(tmp)
        result = work / "profile"
        result.mkdir()
        (work / "input.bam").symlink_to(bam)
        index = next((p for p in (Path(str(bam)+".bai"), bam.with_suffix(".bai")) if p.is_file()), None)
        if index:
            (work / "input.bam.bai").symlink_to(index)
        else:
            subprocess.run(["samtools", "index", "input.bam"], cwd=work, check=True)
        run = lambda args: subprocess.run(args, cwd=work, check=True)
        run(["mhaptools", "convert", "-i", "input.bam", "-c", str(cpg), "-o", "sample.mhap.gz"])
        with gzip.open(work / "sample.mhap.gz", "rt") as fh:
            if not fh.readline().strip():
                raise ValueError("No methylation haplotypes called. Check the hg19 reference, methylation-informative bases and alignment flags.")
        run(["tabix", "-b", "2", "-e", "3", "-p", "bed", "sample.mhap.gz"])
        metrics = ["MM", "PDR", "CHALM", "MHL", "MCR", "MBS", "Entropy"]
        with tempfile.TemporaryFile() as log:
            stat = subprocess.run(["mhaptk", "stat", "--mhapPath", "sample.mhap.gz", "--cpgPath", str(cpg),
                                   "--bedPath", str(bed), "--metrics", *metrics, "--outputFile", str(result/"cgi.tsv")],
                                  cwd=work, stdout=log, stderr=subprocess.STDOUT)
            if stat.returncode:
                log.seek(max(0, log.tell()-4096))
                raise ValueError("mHapTk stat failed\n"+log.read().decode(errors="replace"))
        with (result / "cgi.tsv").open() as fh, bed.open() as intervals:
            table = csv.DictReader(fh, delimiter="\t")
            if not {"chr", "start", "end", "nReads", *metrics}.issubset(table.fieldnames or []):
                raise ValueError("Unrecognized mHapTk statistics header")
            count = 0
            bed_rows = (row for row in csv.reader(intervals, delimiter="\t") if row)
            for region, row in zip_longest(bed_rows, table):
                if (region is None or row is None or len(region)<3 or
                        (region[0], int(region[1]), int(region[2])) !=
                        (row["chr"], int(row["start"]), int(row["end"]))):
                    raise ValueError("mHapTk output does not cover the supplied BED in order")
                count += 1
            if not count:
                raise ValueError("Empty region statistics")
        run(["Rscript", "--vanilla", str(ROOT / "R/motifs.R"), "input.bam", str(result / "motifs.tsv")])
        with (result / "motifs.tsv").open() as source:
            motifs = list(csv.DictReader(source, delimiter="\t"))
        for channel in ("EM", "BPM"):
            valid_total = sum(int(row[channel]) for row in motifs if "N" not in row["kmer"])
            if valid_total == 0:
                raise ValueError("No unambiguous " + channel + " fragments")
            with (result / (channel + ".tsv")).open("w", newline="") as dest:
                writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
                writer.writerow(["kmer", "count", "frequency"])
                for row in motifs:
                    value = int(row[channel])
                    frequency = value / valid_total if "N" not in row["kmer"] else "NA"
                    writer.writerow([row["kmer"], value, frequency])
        for suffix in ("", ".tbi"):
            (work / ("sample.mhap.gz" + suffix)).rename(result / ("sample.mhap.gz" + suffix))
        if out.exists():
            raise ValueError(f"Output already exists: {out}")
        result.rename(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg, help_text in (("bam","Coordinate-sorted PRIMER BAM"), ("cpg","Indexed hg19 CpG annotation (.gz)"),
                           ("bed","Study CpG-island BED"), ("out","New profile directory for one sample")):
        parser.add_argument("--"+arg, required=True, type=Path, help=help_text)
    try:
        preprocess(**vars(parser.parse_args()))
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, str(exc)+"\n")
