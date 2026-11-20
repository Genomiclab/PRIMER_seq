#!/usr/bin/env python3
"""Build cohort matrices from per-sample mHap files and cohort-wide EM/BPM matrices."""
import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import tempfile
from cohorts import COHORTS, ROOT, read_manifest


def prepare(args):
    rows = read_manifest(args.manifest.resolve(strict=True), ("mhap", "cgi"))
    em, bpm = args.em.resolve(strict=True), args.bpm.resolve(strict=True)
    cpg = args.cpg.resolve(strict=True) if args.cpg else None
    if any(row["kind"] == "mhap" for row in rows) and cpg is None:
        raise ValueError("mHap input requires --cpg and its adjacent .tbi index")
    out = args.out.absolute()
    if out.exists() or args.jobs < 1:
        raise ValueError("Supply a new output directory and positive job count")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".processed-", dir=out.parent) as tmp:
        work = Path(tmp) / "run"
        for name in ("profiles", "manifests", "matrices"):
            (work / name).mkdir(parents=True)

        def summarize(row):
            if row["kind"] == "cgi":
                return row["src"]
            dest = work / "profiles" / row["sample_id"]
            subprocess.run([sys.executable, str(ROOT / "scripts/step1_mhap_to_cgi.py"),
                            "--mhap", str(row["src"]), "--cpg", str(cpg), "--bed", str(args.bed),
                            "--out", str(dest)], check=True)
            return dest / "cgi.tsv"

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            paths = dict(zip((r["sample_id"] for r in rows), pool.map(summarize, rows)))
        for cohort in COHORTS:
            subset = [r for r in rows if r["cohort"] == cohort]
            if not subset:
                continue
            columns = ["sample_id", "sex", "cgi"]
            columns += [k for k in ("participant_id", "EM_valid_fraction", "BPM_valid_fraction") if k in subset[0]]
            if cohort == "train":
                columns += ["group"]
                if "methylation_block" in subset[0]:
                    if not all(r["methylation_block"] for r in subset):
                        raise ValueError("Fill methylation_block for every training sample or omit the column")
                    columns += ["methylation_block"]
            manifest = work / "manifests" / (cohort + ".tsv")
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in subset:
                    result = {key: row.get(key, "") for key in columns}
                    path = paths[row["sample_id"]]
                    result["cgi"] = os.path.relpath(path, manifest.parent) if row["kind"] == "mhap" else str(path)
                    writer.writerow(result)
            subprocess.run(["Rscript", "--vanilla", str(ROOT / "scripts/step2_build_matrices.R"),
                            "--cgi", str(manifest), "--em", str(em), "--bpm", str(bpm),
                            "--motif_format", args.motif_format, "--layout", args.layout,
                            "--out", str(work / "matrices" / cohort)], check=True)
        if out.exists():
            raise ValueError("Output already exists")
        work.rename(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("manifest", "em", "bpm", "out"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--cpg", type=Path)
    parser.add_argument("--bed", type=Path, default=ROOT / "config/hg19_cpgIsland.bed")
    parser.add_argument("--motif-format", choices=("counts", "canonical-frequencies", "all-frequencies"), default="counts")
    parser.add_argument("--layout", choices=("auto", "samples-by-features", "features-by-samples"), default="auto")
    parser.add_argument("--jobs", type=int, default=1)
    try:
        prepare(parser.parse_args())
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, str(exc) + "\n")
