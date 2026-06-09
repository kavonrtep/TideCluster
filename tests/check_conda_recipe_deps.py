#!/usr/bin/env python3
"""Assert `conda/tidecluster/meta.yaml`'s `requirements.run` block stays in
sync with `conda-deps.txt` — both that every package is present AND that every
*version pin* matches.

Two classes of drift this guards against:
- **Missing package** (the 1.10.1 packaging slip): `kitehor` was added to
  `conda-deps.txt` but never reached the recipe, so the conda package shipped
  without it.
- **Stale pin** (the 1.14.0 slip): `conda-deps.txt` moved kitehor to `0.13.2`
  and r-igraph to `=2.0.3` while the recipe still said `kitehor=0.12.0` and an
  unpinned `r-igraph`, so `mamba install tidecluster` pulled the wrong KITE
  backend and an unpinned (non-deterministic) igraph.

`conda-deps.txt` is the single source of truth. For every package it pins
(`pkg=X` / `pkg>=X` / …), the recipe MUST carry the identical spec. Packages
listed without a pin only need to be present by name (either side may stay
unpinned). Invoked from `tests/smoke.sh` and any release-gate test.

Implementation notes
- Plain-text parsing (no `pyyaml` dep); the recipe's `requirements.run` is a
  flat list, easy to parse with `re`.
- Comments (`#`) are stripped from both files before parsing.
"""
import re, sys, pathlib

ROOT  = pathlib.Path(__file__).resolve().parent.parent
DEPS  = ROOT / "conda-deps.txt"
META  = ROOT / "conda" / "tidecluster" / "meta.yaml"


def _strip_pin(token: str) -> str:
    """`pkg=1.2.3` / `pkg>=1.2` / `pkg` → `pkg`."""
    return re.split(r"[=<>!~]", token, 1)[0].strip()


def _is_pinned(spec: str) -> bool:
    return _strip_pin(spec) != spec.strip()


def read_conda_deps(path: pathlib.Path) -> "dict[str, str]":
    """name(lower) -> full spec (e.g. 'kitehor' -> 'kitehor=0.13.2')."""
    specs: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        specs[_strip_pin(line).lower()] = line
    return specs


def read_recipe_run(path: pathlib.Path) -> "dict[str, str]":
    """name(lower) -> full spec, from the recipe's requirements.run block."""
    text = path.read_text()
    m = re.search(
        r"^requirements:\s*$\s*^\s*run:\s*$([\s\S]*?)(?=^[A-Za-z])",
        text, re.M)
    if not m:
        sys.exit("FAIL: could not locate `requirements.run` block in "
                 f"{path}")
    specs: dict[str, str] = {}
    for raw in m.group(1).splitlines():
        line = raw.split("#", 1)[0].strip()
        mm = re.match(r"-\s*([A-Za-z0-9_.\-]+(?:\s*[=<>!~][^#]*)?)", line)
        if mm:
            spec = mm.group(1).strip()
            specs[_strip_pin(spec).lower()] = spec
    return specs


def main() -> int:
    declared = read_conda_deps(DEPS)
    in_recipe = read_recipe_run(META)

    missing = [spec for name, spec in declared.items() if name not in in_recipe]
    mismatched = [
        (spec, in_recipe[name])
        for name, spec in declared.items()
        if name in in_recipe and _is_pinned(spec) and in_recipe[name] != spec
    ]

    if missing:
        print("FAIL: packages in conda-deps.txt but missing from "
              "conda/tidecluster/meta.yaml requirements.run:\n  "
              + "\n  ".join(missing))
        print("\nAdd them to the recipe so `mamba install tidecluster` "
              "actually pulls them in.")
        return 1
    if mismatched:
        print("FAIL: version pin drift between conda-deps.txt and "
              "conda/tidecluster/meta.yaml requirements.run:")
        for want, got in mismatched:
            print(f"  conda-deps.txt: {want:<22}  recipe: {got}")
        print("\nconda-deps.txt is the source of truth — update the recipe "
              "pins to match, else `mamba install tidecluster` ships the "
              "wrong version.")
        return 1

    print(f"check_conda_recipe_deps: all {len(declared)} packages present and "
          "every conda-deps.txt pin matches meta.yaml requirements.run ✓")
    return 0


if __name__ == "__main__":
    sys.exit(main())
