"""Generate per-tissue input files from the RatGTEx sample table so that FASTQs
can be downloaded from the SRA by the pipeline.

The RatGTEx sample table
(https://ratgtex.org/data/v4/ref/RatGTEx_samples.v4.tsv) maps each RatGTEx
sample to its SRA accession. Its columns are:

    rat_id  tissue  GEO_accession  BioSample_accession  SRA_accession

`SRA_accession` is generally an SRA *experiment* accession (SRX...), which may
contain one or more sequencing *runs* (SRR...). This script resolves each
experiment to its run accession(s) and library layout (paired- vs single-end)
using the public ENA portal API, then writes, for each requested tissue:

    {version}/{tissue}/rat_ids.txt   one rat ID per line (only rats with runs)
    {version}/{tissue}/fastq_map.txt tab-delimited FASTQ-to-rat map used by the
                                     pipeline. For paired-end runs each row is
                                     `{run}_1.fastq.gz  {run}_2.fastq.gz  rat_id`;
                                     for single-end runs each row is
                                     `{run}.fastq.gz  rat_id`. Paths are relative
                                     to the tissue's `fastq_path`, so the
                                     download rules write them there.
    {version}/{tissue}/sra_runs.tsv  record of rat_id/experiment/run/layout used
                                     to build the FASTQ map (for reference).

The FASTQ files themselves are produced from these run accessions by the
`sra_fastq_paired` / `sra_fastq_single` Snakemake rules (steps/download.smk),
which run `prefetch` + `fasterq-dump` from the SRA Toolkit.

Run this on a machine with internet access before running Snakemake, e.g.:

    python3 scripts/setup/setup_sra.py \
        --samples RatGTEx_samples.v4.tsv \
        --version v4 \
        --tissue Liver
"""

import argparse
import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RUN_PREFIXES = ("SRR", "ERR", "DRR")


def read_samples(path: Path) -> list[dict]:
    """Read the RatGTEx sample table into a list of row dicts."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"rat_id", "tissue", "SRA_accession"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}. "
                f"Found columns: {', '.join(reader.fieldnames or [])}"
            )
        return list(reader)


def resolve_runs(accession: str, retries: int = 3) -> list[tuple[str, bool]]:
    """Resolve an SRA accession to its run accessions and layout via ENA.

    Returns a list of (run_accession, is_paired) tuples. If the accession is
    already a run accession and cannot be looked up, it is returned as-is with
    paired-end assumed unknown (treated as single-end); prefer resolving so the
    layout is correct.
    """
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,library_layout",
            "format": "tsv",
        }
    )
    url = f"{ENA_FILEREPORT}?{query}"
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                text = resp.read().decode()
            break
        except Exception as err:  # noqa: BLE001 - report and retry network errors
            last_err = err
            time.sleep(2 ** attempt)
    else:
        raise RuntimeError(f"Failed to query ENA for {accession}: {last_err}")

    runs = []
    lines = text.splitlines()
    for line in lines[1:]:  # skip header
        fields = line.split("\t")
        if len(fields) < 2 or not fields[0]:
            continue
        run_acc = fields[0].strip()
        is_paired = fields[1].strip().upper() == "PAIRED"
        runs.append((run_acc, is_paired))
    return runs


def write_tissue(version: str, tissue: str, rows: list[dict], no_resolve: bool) -> None:
    """Write rat_ids.txt, fastq_map.txt, and sra_runs.tsv for one tissue."""
    out_dir = Path(version) / tissue
    out_dir.mkdir(parents=True, exist_ok=True)

    rat_ids: list[str] = []
    fastq_rows: list[list[str]] = []
    run_records: list[list[str]] = []

    for row in rows:
        rat_id = row["rat_id"].strip()
        accession = (row.get("SRA_accession") or "").strip()
        if not accession:
            print(f"WARNING: {tissue}/{rat_id} has no SRA_accession; skipping", file=sys.stderr)
            continue

        if no_resolve or accession.startswith(RUN_PREFIXES):
            # Treat the accession itself as a run; layout unknown without lookup.
            runs = [(accession, None)]
        else:
            runs = resolve_runs(accession)
            if not runs:
                print(f"WARNING: no runs found for {tissue}/{rat_id} ({accession}); skipping", file=sys.stderr)
                continue

        rat_ids.append(rat_id)
        for run_acc, is_paired in runs:
            if is_paired is None:
                # Unknown layout: default to single-end and warn.
                print(
                    f"WARNING: layout unknown for {run_acc} ({tissue}/{rat_id}); "
                    "assuming single-end. Resolve via ENA to confirm.",
                    file=sys.stderr,
                )
                is_paired = False
            if is_paired:
                fastq_rows.append([f"{run_acc}_1.fastq.gz", f"{run_acc}_2.fastq.gz", rat_id])
            else:
                fastq_rows.append([f"{run_acc}.fastq.gz", rat_id])
            run_records.append([rat_id, accession, run_acc, "PAIRED" if is_paired else "SINGLE"])

    (out_dir / "rat_ids.txt").write_text("\n".join(rat_ids) + "\n" if rat_ids else "")
    with open(out_dir / "fastq_map.txt", "w") as f:
        for fields in fastq_rows:
            f.write("\t".join(fields) + "\n")
    with open(out_dir / "sra_runs.tsv", "w") as f:
        f.write("rat_id\tSRA_accession\trun_accession\tlayout\n")
        for rec in run_records:
            f.write("\t".join(rec) + "\n")

    print(f"{tissue}: {len(rat_ids)} rats, {len(fastq_rows)} FASTQ (pairs/files) -> {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", required=True, type=Path, help="Path to RatGTEx_samples.v4.tsv")
    parser.add_argument("--version", required=True, help="RatGTEx version, e.g. v4")
    parser.add_argument(
        "--tissue",
        action="append",
        dest="tissues",
        help="Tissue to set up (repeatable). Omit to set up all tissues in the table.",
    )
    parser.add_argument(
        "--no-resolve",
        action="store_true",
        help="Do not query ENA; treat SRA_accession as a run accession as-is (offline).",
    )
    args = parser.parse_args()

    rows = read_samples(args.samples)
    tissues = args.tissues or sorted({row["tissue"] for row in rows})
    by_tissue: dict[str, list[dict]] = {t: [] for t in tissues}
    for row in rows:
        if row["tissue"] in by_tissue:
            by_tissue[row["tissue"]].append(row)

    for tissue in tissues:
        if not by_tissue[tissue]:
            print(f"WARNING: no rows for tissue {tissue} in {args.samples}", file=sys.stderr)
            continue
        write_tissue(args.version, tissue, by_tissue[tissue], args.no_resolve)


if __name__ == "__main__":
    main()
