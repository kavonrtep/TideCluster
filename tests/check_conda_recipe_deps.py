#!/usr/bin/env python3
"""Assert every package listed in `conda-deps.txt` is also present in
`conda/tidecluster/meta.yaml`'s `requirements.run` block.

Catches the kind of drift that produced the 1.10.1 packaging slip:
`kitehor=0.12.0` was added to `conda-deps.txt` when the kitehor
integration landed but never made it into the recipe, so the conda
package shipped without it. Invoked from `tests/smoke.sh` and from any
release-gate test that wants a fast check.

Implementation notes
- Compares only package *names* (ignoring version pins) — version
  pins legitimately differ across files (one might float, the other
  pin to a stable upstream).
- Plain-text YAML parsing to avoid pulling in `pyyaml` as a smoke-test
  dep; the recipe's `requirements.run` is a flat list and easy to
  parse with `re`.
- Comments (`#`) are stripped from both files before name extraction.
"""
import re, sys, pathlib

ROOT  = pathlib.Path(__file__).resolve().parent.parent
DEPS  = ROOT / "conda-deps.txt"
META  = ROOT / "conda" / "tidecluster" / "meta.yaml"


def _strip_pin(token: str) -> str:
    """`pkg=1.2.3` / `pkg>=1.2` / `pkg` → `pkg`."""
    return re.split(r"[=<>!~]", token, 1)[0].strip()


def read_conda_deps(path: pathlib.Path) -> list[str]:
    pkgs: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        pkgs.append(_strip_pin(line))
    return pkgs


def read_recipe_run(path: pathlib.Path) -> set[str]:
    text = path.read_text()
    # Find `requirements:` then its `run:` child block (until the next
    # top-level YAML key — `[a-z]` at column 0 — or end of file).
    m = re.search(
        r"^requirements:\s*$\s*^\s*run:\s*$([\s\S]*?)(?=^[A-Za-z])",
        text, re.M)
    if not m:
        sys.exit("FAIL: could not locate `requirements.run` block in "
                 f"{path}")
    block = m.group(1)
    pkgs: set[str] = set()
    for raw in block.splitlines():
        line = raw.split("#", 1)[0].strip()
        # List item: `- pkgname` or `- pkgname=1.2.3`.
        mm = re.match(r"-\s*([A-Za-z0-9_.\-]+)", line)
        if mm:
            pkgs.add(_strip_pin(mm.group(1)).lower())
    return pkgs


def main() -> int:
    declared = read_conda_deps(DEPS)
    in_recipe = read_recipe_run(META)
    missing = [p for p in declared if p.lower() not in in_recipe]
    if missing:
        print(
            "FAIL: packages in conda-deps.txt but missing from "
            "conda/tidecluster/meta.yaml requirements.run:\n  "
            + "\n  ".join(missing))
        print(
            "\nAdd them to the recipe so `mamba install tidecluster` "
            "actually pulls them in. See conda/tidecluster/meta.yaml.")
        return 1
    print(
        f"check_conda_recipe_deps: all {len(declared)} packages in "
        "conda-deps.txt also present in meta.yaml requirements.run ✓")
    return 0


if __name__ == "__main__":
    sys.exit(main())
