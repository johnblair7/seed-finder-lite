"""
Cluster-friendly utilities: split guides into chunks for parallel jobs, merge results.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List, Union


def split_guides(
    input_path: Union[str, Path],
    n_chunks: int,
    output_dir: Union[str, Path],
    *,
    prefix: str = "chunk",
) -> List[Path]:
    """
    Split a guide CSV into N chunk files for parallel cluster jobs.

    Each chunk is a valid CSV with the same header. Use with job arrays:
    run targetfinder chunk_001.csv ... -o results_001.csv, etc., then merge_results().

    Parameters
    ----------
    input_path : path to CSV/TSV with guide_ID and sequence
    n_chunks : number of chunk files to create
    output_dir : directory to write chunk_*.csv
    prefix : filename prefix (chunk_000.csv, chunk_001.csv, ...)

    Returns
    -------
    List of paths to chunk files.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    delim = "\t" if input_path.suffix in (".tsv", ".txt") else ","
    with open(input_path) as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        rows = list(reader)

    total = len(rows)
    if total == 0:
        raise ValueError("No data rows in input.")
    chunk_size = (total + n_chunks - 1) // n_chunks
    out_paths = []

    for i in range(n_chunks):
        start = i * chunk_size
        end = min(start + chunk_size, total)
        if start >= total:
            break
        chunk_path = output_dir / f"{prefix}_{i:03d}.csv"
        with open(chunk_path, "w", newline="") as out:
            w = csv.writer(out, delimiter=",")
            w.writerow(header)
            w.writerows(rows[start:end])
        out_paths.append(chunk_path)

    return out_paths


def merge_results(
    csv_paths: List[Union[str, Path]],
    output_path: Union[str, Path],
) -> Path:
    """
    Merge multiple off-target result CSVs (e.g. from cluster chunk runs) into one.

    Assumes all have the same header. Writes one header and concatenates data rows.

    Parameters
    ----------
    csv_paths : paths to result CSVs in any order
    output_path : merged output CSV

    Returns
    -------
    Path to the merged file.
    """
    output_path = Path(output_path)
    cols = None
    writer = None

    with open(output_path, "w", newline="") as out:
        for path in csv_paths:
            path = Path(path)
            if not path.exists():
                continue
            with open(path) as f:
                reader = csv.DictReader(f)
                if cols is None:
                    cols = reader.fieldnames
                    writer = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
                    writer.writeheader()
                for row in reader:
                    if writer:
                        writer.writerow(row)

    return output_path


def merge_results_simple(
    csv_paths: List[Union[str, Path]],
    output_path: Union[str, Path],
) -> Path:
    """
    Merge result CSVs: same header, concatenate data rows. Simple implementation.
    """
    output_path = Path(output_path)
    written_header = False

    with open(output_path, "w", newline="") as out:
        for path in csv_paths:
            path = Path(path)
            if not path.exists():
                continue
            with open(path) as f:
                lines = f.readlines()
            if not lines:
                continue
            if not written_header:
                out.writelines(lines)
                written_header = True
            else:
                # skip header line
                out.writelines(lines[1:])

    return output_path
