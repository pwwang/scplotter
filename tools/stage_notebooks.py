#!/usr/bin/env python3
"""Stage the Giotto notebooks at the repo root so they can be executed by
`jupyter nbconvert --execute` from the repo root against an external data
directory, without editing the committed notebooks.

In every code cell, the hard-coded conda python path is replaced with the
value of the GIOTTO_PYTHON environment variable.

usage: python3 tools/stage_notebooks.py [--dry-run]
"""

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC_DIR = ROOT / "notebooks" / "spatial"

OLD = ('python_path <- file.path(Sys.getenv("HOME"), "miniconda3", '
       '"envs", "giotto_env", "bin", "python")')
NEW = 'python_path <- Sys.getenv("GIOTTO_PYTHON")'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="print what would be staged without writing")
    args = parser.parse_args()

    staged = 0
    rewrote = 0
    for src in sorted(SRC_DIR.glob("Giotto.*.ipynb")):
        nb = json.loads(src.read_text())
        k = 0
        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            for i, line in enumerate(cell["source"]):
                if line.rstrip("\r\n") == OLD:
                    cell["source"][i] = NEW + ("\n" if line.endswith("\n") else "")
                    k += 1
        dst = ROOT / src.name
        if args.dry_run:
            print(f"{src.relative_to(ROOT)} -> {dst.relative_to(ROOT)} (rewrote {k} cells)")
        else:
            dst.write_text(json.dumps(nb, indent=1))
            print(f"{src.relative_to(ROOT)} -> {dst.relative_to(ROOT)} (rewrote {k} cells)")
        staged += 1
        rewrote += k
    print(f"staged {staged} notebooks, rewrote {rewrote} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
