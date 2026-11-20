#!/usr/bin/env python3
"""Calculate CpG-island methylation features from an existing mHap file."""
import argparse
import csv
import gzip
from itertools import zip_longest
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("MM", "PDR", "CHALM", "MHL", "MCR", "MBS", "Entropy")


def summarize(mhap, cpg, bed, out):
    mhap, cpg, bed = [p.resolve(strict=True) for p in (mhap, cpg, bed)]
    out = out.absolute()
    if out.exists():
        raise ValueError("Use a new output directory: " + str(out))
    if not Path(str(cpg) + ".tbi").is_file():
        raise ValueError("The CpG reference requires an adjacent .tbi index")
    for tool in ("mhaptk", "tabix"):
        if not shutil.which(tool):
            raise ValueError(tool + " must be on PATH")
    with gzip.open(mhap, "rt") as handle:
        row = handle.readline().strip().split("\t")
        if len(row) != 6:
            raise ValueError("Expected a nonempty six-column, bgzip-compressed mHap file")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".mhap-", dir=out.parent) as tmp:
        work = Path(tmp)
        (work / "sample.mhap.gz").symlink_to(mhap)
        index = Path(str(mhap) + ".tbi")
        if index.is_file():
            (work / "sample.mhap.gz.tbi").symlink_to(index)
        else:
            subprocess.run(["tabix", "-b", "2", "-e", "3", "-p", "bed", "sample.mhap.gz"], cwd=work, check=True)
        profile = work / "profile"
        profile.mkdir()
        with tempfile.TemporaryFile() as log:
            result = subprocess.run(["mhaptk", "stat", "--mhapPath", "sample.mhap.gz", "--cpgPath", str(cpg),
                                     "--bedPath", str(bed), "--metrics", *METRICS,
                                     "--outputFile", str(profile / "cgi.tsv")],
                                    cwd=work, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                log.seek(max(0, log.tell() - 4096))
                raise ValueError("mHapTk stat failed\n" + log.read().decode(errors="replace"))
        with (profile / "cgi.tsv").open() as handle, bed.open() as intervals:
            table = csv.DictReader(handle, delimiter="\t")
            if not {"chr", "start", "end", "nReads", *METRICS}.issubset(table.fieldnames or []):
                raise ValueError("Unexpected mHapTk statistics header")
            n = 0
            for region, row in zip_longest((r for r in csv.reader(intervals, delimiter="\t") if r), table):
                if region is None or row is None or tuple(region[:3]) != (row["chr"], row["start"], row["end"]):
                    raise ValueError("Statistics do not cover the supplied BED in order")
                n += 1
            if not n:
                raise ValueError("Empty region statistics")
        if out.exists():
            raise ValueError("Output already exists")
        profile.rename(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhap", required=True, type=Path)
    parser.add_argument("--cpg", required=True, type=Path)
    parser.add_argument("--bed", type=Path, default=ROOT / "config/hg19_cpgIsland.bed")
    parser.add_argument("--out", required=True, type=Path)
    try:
        summarize(**vars(parser.parse_args()))
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, str(exc) + "\n")
