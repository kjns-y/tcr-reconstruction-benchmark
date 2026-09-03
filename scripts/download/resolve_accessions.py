#!/usr/bin/env python3
"""Resolve GEO GSM accessions to SRA runs using current NCBI metadata."""

from __future__ import annotations

import argparse
import csv
import re
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full"
RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={experiment}"


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def fetch_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "tcrbench/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Unable to retrieve NCBI metadata: {url}: {exc}") from exc


def resolve_gsm(gsm: str) -> List[Dict[str, str]]:
    soft = fetch_text(GEO_URL.format(gsm=gsm))
    title_match = re.search(r"^!Sample_title = (.+)$", soft, flags=re.MULTILINE)
    srx_match = re.search(r"^!Sample_relation = SRA: .*term=(SRX\d+)$", soft, flags=re.MULTILINE)
    if not srx_match:
        raise RuntimeError(f"No SRA experiment relation found for {gsm}; refusing to guess an SRR.")
    experiment = srx_match.group(1)
    title = title_match.group(1).strip() if title_match else ""
    runinfo = fetch_text(RUNINFO_URL.format(experiment=experiment))
    rows = list(csv.DictReader(runinfo.splitlines()))
    if not rows or not all(row.get("Run", "").startswith("SRR") for row in rows):
        raise RuntimeError(f"No SRR runs returned by NCBI for {gsm} / {experiment}.")
    output: List[Dict[str, str]] = []
    for row in rows:
        output.append(
            {
                "gsm": gsm,
                "title": title,
                "experiment": experiment,
                "run": row.get("Run", ""),
                "spots": row.get("spots", ""),
                "bases": row.get("bases", ""),
                "size_MB": row.get("size_MB", ""),
                "library_layout": row.get("LibraryLayout", ""),
                "platform": row.get("Platform", ""),
                "model": row.get("Model", ""),
                "download_path": row.get("download_path", ""),
            }
        )
    return output


def write_rows(path: Path, rows: Iterable[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "gsm", "title", "experiment", "run", "spots", "bases", "size_MB",
        "library_layout", "platform", "model", "download_path",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gsm", nargs="+", default=["GSM3148575", "GSM3148580"])
    parser.add_argument("--output-dir", type=Path, default=project_root() / "data/metadata")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for gsm in args.gsm:
        rows = resolve_gsm(gsm)
        output = args.output_dir / f"{gsm}_accessions.tsv"
        write_rows(output, rows)
        print(f"{gsm}: {len(rows)} run(s) -> {output}")


if __name__ == "__main__":
    main()
