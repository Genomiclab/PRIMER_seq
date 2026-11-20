#!/usr/bin/env python3
"""Run the numbered workflow on the supplied synthetic BAM and three-cancer cohort."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def demo(out, threads):
    out = out.absolute()
    if out.exists() or threads < 1:
        raise ValueError("Supply a new output directory and positive thread count")
    out.mkdir(parents=True)

    def run(label, command):
        print("\n" + label, flush=True)
        subprocess.run([str(x) for x in command], cwd=ROOT, check=True)

    def r(step, *args):
        return ["Rscript", "--vanilla", ROOT / "scripts" / step, *args]

    run("Checking dependencies", r("check_environment.R", "all"))
    bam = ROOT / "examples" / "bam"
    run("Step 1: BAM -> mHap, five methylation features, EM and BPM",
        [sys.executable, ROOT / "scripts/step1_bam_to_features.py", "--bam", bam / "synthetic.bam",
         "--cpg", bam / "hg19_demo_CpG.gz", "--bed", ROOT / "config/hg19_cpgIsland.bed", "--out", out / "step1"])
    manifest = out / "bam_sample.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "sex", "profile"])
        writer.writerow(["SIM_BAM_001", "U", "step1"])
    run("Step 2: assemble the BAM-derived feature matrices",
        r("step2_build_matrices.R", "--profiles", manifest, "--out", out / "matrices/bam"))
    tables = out / "tables"
    run("Prepare synthetic matrices for the complete marker panel",
        ["Rscript", "--vanilla", ROOT / "examples/make_study_demo.R", "--out", tables])
    for cohort in ("train", "validation", "blind"):
        run("Step 2: import the synthetic " + cohort + " cohort",
            r("step2_build_matrices.R", "--tables", tables / cohort, "--out", out / "matrices" / cohort))
    run("Step 3: train all five learners and both second-layer stacks",
        r("step3_train_ensemble.R", "--train", out / "matrices/train",
          "--out", out / "model", "--preset", "demo", "--threads", threads))
    for cohort in ("bam", "validation", "blind"):
        run("Step 4: predict " + cohort + " with the newly trained model",
            r("step4_predict.R", "--data", out / "matrices" / cohort, "--model", out / "model",
              "--out", out / "predictions" / cohort, "--cohort", "blind" if cohort == "blind" else "validation"))
        if cohort != "bam":
            run("Step 5: join withheld synthetic labels and evaluate " + cohort,
                r("step5_evaluate.R", "--predictions", out / "predictions" / cohort / "predictions.tsv",
                  "--truth", tables / cohort / "truth.tsv", "--out", out / "reports" / cohort))
    (out / "COMPLETED.txt").write_text("Synthetic demonstration completed. These results are not clinical performance estimates.\n")
    print("\nCompleted: " + str(out), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "work/demo")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    try:
        demo(args.out, args.threads)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, str(exc) + "\n")
