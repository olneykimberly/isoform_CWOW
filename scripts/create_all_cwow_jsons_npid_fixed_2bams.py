#!/usr/bin/env python3
"""
Generate one GrandCanyon JSON per CWOW sample by calling fastq_to_json.py.
The GrandCanyon pipeline is part of the available TGen in-house workflows: https://github.com/tgen

This wrapper is specific to my CWOW long read project. It contains hard-coded paths and would need to be modified if used for a different project.

For each cohort sample, this wrapper:
  1. Reads sex from the cohort CSV. Sex was already confirmed from both the short-read data and SNP array data.
  2. Discovers that sample's PacBio Kinnex *.hifi.u.bam files.
  3. Requires the expected number of BAMs/flowcells (default: 2).
  4. Selects the XY male- or XX female-specific reference FASTA.
  5. Adds an external SJ.out.tab entry when one is found.
  6. Calls fastq_to_json.py once for that sample.

Example
-------
python3 create_all_cwow_jsons.py \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \ # full cohort metdata
  --fastq-to-json /path/to/fastq_to_json.py \
  --template /path/to/grandcanyon_template.json \ # see git repo: https://github.com/tgen/GrandCanyon
  --bam-base /tgen_labs/jfryer/cores/tgen \
  --output-dir /path/to/output/jsons \
  --results-base /path/to/store/grandcanyon/results \
  --male-reference /path/to/GRCh38_male.fa \
  --female-reference /path/to/GRCh38_female.fa \
  --gtf /tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation.gtf \
  --sj-root /path/to/STAR_SJ \ # path to short read SJ
  --account tgen-206000 \ # my TGen HPC account information
  --overwrite \ # to overwrite the sample.json if one already exists
  --dry-run # remove and add --submit if wanting to run
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


BAM_SUFFIX = ".hifi.u.bam"
SEX_MAP = {
    "M": "Male",
    "MALE": "Male",
    "F": "Female",
    "FEMALE": "Female",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create one GrandCanyon JSON per CWOW sample."
    )
    parser.add_argument("--cohort-csv", required=True, type=Path)
    parser.add_argument("--fastq-to-json", required=True, type=Path)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument(
        "--bam-base",
        default=Path("/tgen_labs/jfryer/cores/tgen"),
        type=Path,
        help="Root containing PacBio run directories.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--results-base",
        required=True,
        type=Path,
        help="Base archive directory; each sample gets its own subdirectory.",
    )
    parser.add_argument("--male-reference", required=True, type=Path)
    parser.add_argument("--female-reference", required=True, type=Path)
    parser.add_argument(
        "--gtf",
        required=True,
        type=Path,
        help="Reference GTF annotation used by GrandCanyon RNA processing.",
    )
    parser.add_argument(
        "--sj-root",
        type=Path,
        default=None,
        help="Root containing external short-read SJ.out.tab files.",
    )
    parser.add_argument(
        "--sj-pattern",
        default="{sample}/SJ.out.tab",
        help=(
            "Path beneath --sj-root. Available placeholders: {sample}, "
            "{subject}, and {cohort_id}."
        ),
    )
    parser.add_argument(
        "--sample-column",
        default="Study-Subject-ID",
        help="Cohort CSV column containing IDs such as NA13093-CONTROL.",
    )
    parser.add_argument(
        "--npid-column",
        default="NPID",
        help="Cohort CSV column used to match external SJ.out.tab filenames.",
    )
    parser.add_argument(
        "--sex-column",
        default="Gender",
        help="Cohort CSV column containing Male/Female or M/F.",
    )
    parser.add_argument(
        "--fraction-column",
        default="Specimen_Fraction",
        help="Optional CSV M/F column used to validate --sex-column.",
    )
    parser.add_argument(
        "--expected-bams",
        type=int,
        default=2,
        help="Expected number of uBAMs per sample; use 0 to disable.",
    )
    parser.add_argument(
        "--pipeline",
        default="grandcanyon@development",
    )
    parser.add_argument("--study", default="CWOW")
    parser.add_argument("--assay", default="PKFLR")
    parser.add_argument("--glprep", default="KinnexFullLength")
    parser.add_argument("--gltype", default="RNA")
    parser.add_argument("--platform", default="PACBIO")
    parser.add_argument("--model", default="REVIO")
    parser.add_argument(
        "--account",
        default=None,
        help="Optional SLURM account; otherwise fastq_to_json.py detects it.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Submit each generated JSON to Jetstream Centro.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running fastq_to_json.py.",
    )
    parser.add_argument(
        "--allow-missing-sj",
        action="store_true",
        help="Generate JSON even when an SJ.out.tab file is not found.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate JSON files that already exist.",
    )
    return parser.parse_args()


def normalize_sex(value: object) -> str:
    if pd.isna(value):
        return ""
    return SEX_MAP.get(str(value).strip().upper(), "")


def normalize_cohort_id(value: object) -> str:
    """Convert NA13093_CONTROL or CWOW_NA13093-CONTROL to NA13093-CONTROL."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    text = re.sub(r"^CWOW[_-]", "", text, flags=re.IGNORECASE)
    text = text.replace("_", "-")
    text = re.sub(r"-+", "-", text).strip("-")
    return text.upper()


def project_name(cohort_id: str) -> str:
    """NA13093-CONTROL -> CWOW_NA13093."""
    subject = cohort_id.split("-", 1)[0]
    return f"CWOW_{subject}"


def load_cohort(args: argparse.Namespace) -> dict[str, dict[str, str]]:
    """Load cohort ID, sex, and NPID for every sample."""
    df = pd.read_csv(args.cohort_csv, dtype=str)
    df.columns = df.columns.str.strip()

    required = [args.sample_column, args.sex_column, args.npid_column]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise KeyError(
            f"Missing required CSV column(s): {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    cohort: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        cohort_id = normalize_cohort_id(row[args.sample_column])
        if not cohort_id:
            continue

        sex = normalize_sex(row[args.sex_column])
        if not sex:
            raise ValueError(
                f"Invalid sex value for {cohort_id}: {row[args.sex_column]!r}"
            )

        npid_value = row[args.npid_column]
        if pd.isna(npid_value) or not str(npid_value).strip():
            raise ValueError(f"Missing NPID for {cohort_id}")
        npid = str(npid_value).strip()

        if args.fraction_column in df.columns:
            fraction_sex = normalize_sex(row[args.fraction_column])
            if fraction_sex and fraction_sex != sex:
                raise ValueError(
                    f"Sex conflict for {cohort_id}: "
                    f"{args.sex_column}={sex}, "
                    f"{args.fraction_column}={fraction_sex}"
                )

        previous = cohort.get(cohort_id)
        current = {"sex": sex, "npid": npid}
        if previous and previous != current:
            raise ValueError(
                f"Duplicate cohort ID has conflicting metadata: {cohort_id}"
            )
        cohort[cohort_id] = current

    if not cohort:
        raise RuntimeError("No valid cohort samples were loaded from the CSV.")
    return cohort


def discover_bams(bam_base: Path) -> dict[str, list[Path]]:
    """
    Discover BAMs and group them by full cohort ID.

    Expected sample directory:
      CWOW_NA13093-CONTROL

    Expected BAM:
      CWOW_NA13093-CONTROL_...hifi.u.bam
    """
    grouped: dict[str, list[Path]] = {}

    for bam in bam_base.rglob(f"*{BAM_SUFFIX}"):
        if not bam.is_file():
            continue

        folder = bam.parent.name
        match = re.fullmatch(r"CWOW_(NA\d+-.+)", folder, flags=re.IGNORECASE)
        if not match:
            continue

        cohort_id = normalize_cohort_id(match.group(1))
        if not bam.name.upper().startswith(f"CWOW_{cohort_id}_"):
            continue

        grouped.setdefault(cohort_id, []).append(bam.resolve())

    for cohort_id in grouped:
        grouped[cohort_id] = sorted(set(grouped[cohort_id]))

    return grouped


def resolve_sj_tab(
    args: argparse.Namespace,
    npid: str,
    sex: str,
) -> Path | None:
    """
    Resolve an SJ.out.tab using the cohort CSV NPID.

    Expected examples:
      NA20-446_STAR_XX.bamSJ.out.tab
      NA19-060_STAR_XY.bamSJ.out.tab
    """
    if args.sj_root is None:
        return None

    sex_code = "XY" if sex == "Male" else "XX"
    exact = args.sj_root / f"{npid}_STAR_{sex_code}.bamSJ.out.tab"
    if exact.is_file():
        return exact.resolve()

    # Fallback: tolerate an unexpected or missing XX/XY designation, but only
    # accept the result when exactly one file matches this NPID.
    matches = sorted(args.sj_root.glob(f"{npid}_STAR_*.bamSJ.out.tab"))
    if len(matches) == 1:
        return matches[0].resolve()
    if len(matches) > 1:
        print(
            f"WARNING: multiple SJ files found for NPID {npid}: "
            + ", ".join(str(path) for path in matches),
            file=sys.stderr,
        )
    return None


def input_prefixes(bams: list[Path]) -> list[str]:
    """
    fastq_to_json.py concatenates each input prefix with '*.bam'.
    Use an exact prefix through the stable beginning of each BAM filename.
    """
    return [str(bam.parent / "CWOW") for bam in bams]


def build_command(
    args: argparse.Namespace,
    cohort_id: str,
    sex: str,
    bams: list[Path],
    sj_tab: Path | None,
) -> list[str]:
    sample = project_name(cohort_id)
    reference = (
        args.male_reference if sex == "Male" else args.female_reference
    )

    constants = {
        "grandcanyon": {
            "grch38": {
                "reference_fasta": str(reference.resolve()),
                "reference_fasta_no_alt": str(reference.resolve()),
                "gtf": str(args.gtf.resolve()),
            }
        }
    }

    cmd = [
        sys.executable,
        str(args.fastq_to_json.resolve()),
        "-n", sample,
        "-P", str((args.results_base / sample).resolve()),
        "-p", args.pipeline,
        "-t", str(args.template.resolve()),
        "-i", *input_prefixes(bams),
        "--ubam",
        "-a", args.assay,
        "-gp", args.glprep,
        "-gt", args.gltype,
        "--platform", args.platform,
        "--model", args.model,
        "-s", args.study,
        "--resplit", "_",
        "-lbf", "8",
        "-snf", "1-2",
        "-smf", "1-2",
        "-idf", "9,11,10,8",
        "-o", str(args.output_dir.resolve()),
        "--sex", sex,
        "--constants", json.dumps(constants),
    ]

    if sj_tab is not None:
        external = {
            sample: {
                "sj_tab": str(sj_tab)
            }
        }
        cmd.extend(["--external", json.dumps(external)])

    if args.account:
        cmd.extend(["-A", args.account])
    if args.submit:
        cmd.append("--submit")

    return cmd


def shell_join(cmd: list[str]) -> str:
    import shlex
    return " ".join(shlex.quote(part) for part in cmd)


def validate_paths(args: argparse.Namespace) -> None:
    required_files = [
        args.cohort_csv,
        args.fastq_to_json,
        args.template,
        args.male_reference,
        args.female_reference,
        args.gtf,
    ]
    missing_files = [path for path in required_files if not path.is_file()]
    if missing_files:
        formatted = "\n  ".join(str(path) for path in missing_files)
        raise FileNotFoundError(f"Required file(s) not found:\n  {formatted}")

    if not args.bam_base.is_dir():
        raise NotADirectoryError(f"BAM base does not exist: {args.bam_base}")

    if args.sj_root is not None and not args.sj_root.is_dir():
        raise NotADirectoryError(f"SJ root does not exist: {args.sj_root}")


def main() -> None:
    args = parse_args()
    validate_paths(args)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.results_base.mkdir(parents=True, exist_ok=True)

    cohort = load_cohort(args)
    discovered = discover_bams(args.bam_base)

    generated = 0
    skipped_existing = 0
    missing_bams: list[str] = []
    wrong_bam_count: list[tuple[str, int]] = []
    missing_sj: list[str] = []
    failures: list[tuple[str, int]] = []

    for cohort_id, metadata in sorted(cohort.items()):
        sex = metadata["sex"]
        npid = metadata["npid"]
        sample = project_name(cohort_id)
        bams = discovered.get(cohort_id, [])

        if not bams:
            missing_bams.append(cohort_id)
            continue

        if args.expected_bams and len(bams) != args.expected_bams:
            wrong_bam_count.append((cohort_id, len(bams)))
            continue

        output_json = args.output_dir / f"{sample}.json"
        if output_json.exists() and not args.overwrite:
            print(f"SKIP existing: {output_json}")
            skipped_existing += 1
            continue

        sj_tab = resolve_sj_tab(args, npid, sex)
        if args.sj_root is not None and sj_tab is None:
            missing_sj.append(f"{cohort_id} (NPID={npid})")
            # Missing SJ data is a warning; still generate the long-read JSON.

        cmd = build_command(args, cohort_id, sex, bams, sj_tab)
        print(
            f"\n[{sample}] cohort={cohort_id}, NPID={npid}, "
            f"sex={sex}, BAMs={len(bams)}, SJ={sj_tab or 'not included'}"
        )
        print(shell_join(cmd))

        if args.dry_run:
            generated += 1
            continue

        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            failures.append((cohort_id, result.returncode))
        else:
            generated += 1

    print("\nSummary")
    print(f"  Cohort samples:          {len(cohort)}")
    print(f"  JSONs generated/planned:{generated:>6}")
    print(f"  Existing JSONs skipped: {skipped_existing}")
    print(f"  Samples without BAMs:   {len(missing_bams)}")
    print(f"  Unexpected BAM counts:  {len(wrong_bam_count)}")
    print(f"  Samples without SJ.tab: {len(missing_sj)}")
    print(f"  Command failures:       {len(failures)}")

    if missing_bams:
        print("\nSamples without BAMs:")
        print("  " + "\n  ".join(missing_bams))

    if wrong_bam_count:
        print("\nSamples with unexpected BAM counts:")
        print("  " + "\n  ".join(
            f"{sample}: {count}" for sample, count in wrong_bam_count
        ))

    if missing_sj:
        print("\nSamples without matching SJ.out.tab:")
        print("  " + "\n  ".join(missing_sj))

    if failures:
        print("\nFailed commands:")
        print("  " + "\n  ".join(
            f"{sample}: return code {code}" for sample, code in failures
        ))
        sys.exit(1)


if __name__ == "__main__":
    main()
