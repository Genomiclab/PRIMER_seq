#!/usr/bin/env python3
"""Combine the detection and tissue-of-origin CSV marker lists without reordering features."""
import argparse
import csv
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ("MM", "PDR", "CHALM", "MHL", "Entropy", "EM", "BPM")
TASKS = (("detection", "cancer_detection"), ("origin", "tissue_of_origin"))


def collect(source, bed):
    with bed.open() as handle:
        regions = {".".join(row[:3]) for row in csv.reader(handle, delimiter="\t") if row}
    merged, counts = [], []
    for task, name in TASKS:
        for channel in CHANNELS:
            motif = channel in ("EM", "BPM")
            filename = f"{channel}_{name}_{'motif' if motif else 'methylation'}_markers.csv"
            field = "motif" if motif else "hg19_CGI_genomic_region"
            with (source / filename).open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                if not {"metric", field}.issubset(reader.fieldnames or []):
                    raise ValueError("Missing columns in " + filename)
                features = []
                for row in reader:
                    value = row[field]
                    if row["metric"] != channel or not value or value != value.strip():
                        raise ValueError("Invalid marker in " + filename)
                    if motif:
                        valid = re.fullmatch(r"[ACGT]{6}", value)
                    else:
                        valid = value in regions and re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])[.]\d+[.]\d+", value)
                    if not valid:
                        raise ValueError("Invalid motif or hg19 autosomal region: " + value)
                    features.append(value)
            if len(features) < 2 or len(features) != len(set(features)):
                raise ValueError("Too few or duplicate markers in " + filename)
            merged.extend((task, channel, value) for value in features)
            counts.append((task, channel, len(features), filename))
    return merged, counts


def build(source, bed, out):
    if out.exists():
        raise ValueError("Use a new output directory")
    merged, counts = collect(source, bed)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".markers-", dir=out.parent) as tmp:
        work = Path(tmp) / "panel"
        work.mkdir()
        for filename, header, rows in (
            ("markers.tsv", ("task", "channel", "feature"), merged),
            ("counts.tsv", ("task", "channel", "markers", "source"), counts),
        ):
            with (work / filename).open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)
        if out.exists():
            raise ValueError("Output already exists")
        work.rename(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "marker_file")
    parser.add_argument("--bed", type=Path, default=ROOT / "config/hg19_cpgIsland.bed")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.input, args.bed, args.out)
    except (ValueError, OSError) as exc:
        parser.exit(1, str(exc) + "\n")
