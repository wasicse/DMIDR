#!/usr/bin/env python3
"""Split a multi-sequence FASTA into one file per sequence."""
from __future__ import annotations

import argparse
from pathlib import Path


def split_fasta(input_path: Path, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    header: str | None = None
    lines: list[str] = []

    def flush() -> None:
        if header is None:
            return
        name = header.lstrip(">").split()[0]
        safe_name = name.replace("/", "_").replace("\\", "_").replace(" ", "_")
        out = output_dir / f"{safe_name}.fasta"
        with out.open("w", encoding="utf-8") as fh:
            fh.write(f">{header.lstrip('>')}\n")
            fh.write("".join(lines))
        written.append(out)

    for line in input_path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(">"):
            flush()
            header = line.lstrip(">").rstrip()
            lines = []
        else:
            lines.append(line)
    flush()
    return written


def main() -> None:
    p = argparse.ArgumentParser(description="Split a multi-FASTA into per-sequence files.")
    p.add_argument("--input", required=True, help="Input multi-FASTA file")
    p.add_argument("--output-dir", required=True, help="Directory to write individual FASTA files")
    args = p.parse_args()

    written = split_fasta(Path(args.input), Path(args.output_dir))
    for path in written:
        print(f"  Wrote {path}")
    print(f"Split {len(written)} sequence(s) to {args.output_dir}")


if __name__ == "__main__":
    main()
