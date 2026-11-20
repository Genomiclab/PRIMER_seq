#!/usr/bin/env python3
"""Prepare cohort matrices from a manifest of BAMs or existing PRIMER profiles."""
import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("train", "validation", "blind")
GROUPS = {"Normal", "Lung", "Colorectum", "Stomach", "Liver", "Thyroid", "Breast",
          "Lymphoma", "Esophagus", "Pancreas", "Bladder", "Prostate", "Ovary"}


def read_manifest(path, input_fields=("bam", "profile")):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not {"sample_id", "cohort", "sex"}.issubset(fields):
            raise ValueError("Expected unique columns including sample_id, cohort and sex")
        rows = list(reader)
    if not rows:
        raise ValueError("The manifest is empty")
    samples, participants, sources = set(), set(), set()
    for line, row in enumerate(rows, 2):
        sid = row["sample_id"]
        if not sid or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid) or sid in samples:
            raise ValueError(f"Line {line}: invalid or duplicate sample_id")
        samples.add(sid)
        participant = row.get("participant_id") if "participant_id" in fields else sid
        if not participant or participant in participants:
            raise ValueError(f"Line {line}: supply one specimen per participant across the study cohorts")
        participants.add(participant)
        if row["sex"] in ("", "NA"):
            row["sex"] = "U"
        if row["cohort"] not in COHORTS or row["sex"] not in ("F", "M", "U"):
            raise ValueError(f"Line {line}: unknown cohort or sex")
        if row["cohort"] == "train" and row.get("group") not in GROUPS:
            raise ValueError(f"Line {line}: training requires a study diagnosis in group")
        kinds = [k for k in input_fields if row.get(k)]
        if len(kinds) != 1:
            raise ValueError(f"Line {line}: specify exactly one of " + ", ".join(input_fields))
        kind = kinds[0]
        src = (path.parent / row[kind]).resolve(strict=True)
        if src in sources or (kind != "profile" and not src.is_file()):
            raise ValueError(f"Line {line}: duplicate or invalid input")
        if kind == "profile" and not all((src/f).is_file() for f in ("cgi.tsv", "motifs.tsv")):
            raise ValueError(f"Line {line}: missing profile files")
        sources.add(src)
        row["kind"], row["src"] = kind, src
    return rows


def prepare(args):
    rows = read_manifest(args.manifest.resolve(strict=True))
    out = args.out.absolute()
    if out.exists() or args.jobs < 1:
        raise ValueError("Supply a new output directory and a positive job count")
    if any(r["kind"] == "bam" for r in rows) and args.cpg is None:
        raise ValueError("BAM processing requires --cpg and its adjacent .tbi index")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".cohorts-", dir=out.parent) as tmp:
        work = Path(tmp) / "run"
        for name in ("profiles", "manifests", "matrices"):
            (work/name).mkdir(parents=True)

        def prepare_profile(row):
            if row["kind"] == "profile":
                return row["src"]
            dest = work / "profiles" / row["sample_id"]
            subprocess.run([sys.executable, str(ROOT/"scripts/step1_bam_to_features.py"), "--bam", str(row["src"]),
                            "--cpg", str(args.cpg), "--bed", str(args.bed), "--out", str(dest)], check=True)
            return dest

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            paths = dict(zip((r["sample_id"] for r in rows), pool.map(prepare_profile, rows)))
        for cohort in COHORTS:
            subset = [r for r in rows if r["cohort"] == cohort]
            if not subset:
                continue
            manifest = work / "manifests" / (cohort+".tsv")
            columns = ["sample_id", "sex", "profile"]
            if "participant_id" in subset[0]:
                columns += ["participant_id"]
            if cohort == "train":
                columns += ["group"]
                if "methylation_block" in subset[0]:
                    if not all(r["methylation_block"] for r in subset):
                        raise ValueError("Fill methylation_block for every training sample or omit the column")
                    columns += ["methylation_block"]
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for row in subset:
                    entry = {k: row.get(k, "") for k in columns}
                    path = paths[row["sample_id"]]
                    entry["profile"] = os.path.relpath(path, manifest.parent) if row["kind"] == "bam" else str(path)
                    writer.writerow(entry)
            subprocess.run(["Rscript", "--vanilla", str(ROOT/"scripts/step2_build_matrices.R"), "--profiles",
                            str(manifest), "--out", str(work/"matrices"/cohort)], check=True)
        if out.exists():
            raise ValueError(f"Output already exists: {out}")
        work.rename(out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True, help="Tab-separated study manifest")
    p.add_argument("--cpg", type=Path, help="Indexed hg19 CpG reference")
    p.add_argument("--bed", type=Path, default=ROOT/"config/hg19_cpgIsland.bed")
    p.add_argument("--out", type=Path, required=True, help="New cohort output directory")
    p.add_argument("--jobs", type=int, default=1, help="Number of BAMs to process in parallel")
    try:
        prepare(p.parse_args())
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        p.exit(1, str(exc)+"\n")
