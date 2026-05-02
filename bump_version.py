"""Bump the qPCR Analyzer version across the codebase.

Single utility for keeping every version literal in sync. The canonical
version is read from ``pyproject.toml``; running this script updates every
known location at once.

Usage::

    python bump_version.py                    # show current version
    python bump_version.py 2.0.0              # set explicit version
    python bump_version.py --bump major       # 1.2.3 -> 2.0.0
    python bump_version.py --bump minor       # 1.2.3 -> 1.3.0
    python bump_version.py --bump patch       # 1.2.3 -> 1.2.4
    python bump_version.py 2.0.0 --check      # dry-run; print diff, don't write

Locations updated:

    pyproject.toml                       version = "X.Y.Z"
    src/qpcr_analyzer/__init__.py        __version__ = "X.Y.Z"
    README.md                            @vX.Y.Z (install-URL example)
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SEMVER = r"\d+\.\d+\.\d+"
SEMVER_FULL = re.compile(rf"^{SEMVER}$")


@dataclass(frozen=True)
class Location:
    path: Path
    pattern: re.Pattern[str]
    template: str  # uses {version} placeholder
    label: str


LOCATIONS: list[Location] = [
    Location(
        path=ROOT / "pyproject.toml",
        pattern=re.compile(rf'(?m)^version\s*=\s*"{SEMVER}"\s*$'),
        template='version = "{version}"',
        label="pyproject.toml :: [project] version",
    ),
    Location(
        path=ROOT / "src" / "qpcr_analyzer" / "__init__.py",
        pattern=re.compile(rf'__version__\s*=\s*"{SEMVER}"'),
        template='__version__ = "{version}"',
        label="src/qpcr_analyzer/__init__.py :: __version__",
    ),
    Location(
        path=ROOT / "README.md",
        pattern=re.compile(rf"pip install qpcr-analyzer=={SEMVER}"),
        template="pip install qpcr-analyzer=={version}",
        label="README.md :: PyPI pinned-version example",
    ),
    Location(
        path=ROOT / "README.md",
        pattern=re.compile(rf"py_qpcr_analyzer\.git@v{SEMVER}"),
        template="py_qpcr_analyzer.git@v{version}",
        label="README.md :: source-install tag example",
    ),
]


def read_current_version() -> str:
    """Read the canonical version from pyproject.toml."""
    src = LOCATIONS[0]
    text = src.path.read_text(encoding="utf-8")
    m = src.pattern.search(text)
    if not m:
        raise SystemExit(f"could not locate version in {src.path}")
    inner = re.search(SEMVER, m.group(0))
    if not inner:
        raise SystemExit(f"version match malformed in {src.path}")
    return inner.group(0)


def bump(current: str, kind: str) -> str:
    parts = [int(x) for x in current.split(".")]
    if kind == "major":
        parts = [parts[0] + 1, 0, 0]
    elif kind == "minor":
        parts = [parts[0], parts[1] + 1, 0]
    elif kind == "patch":
        parts = [parts[0], parts[1], parts[2] + 1]
    else:
        raise SystemExit(f"unknown bump kind: {kind}")
    return ".".join(str(p) for p in parts)


def update_location(loc: Location, new_version: str, *, write: bool) -> int:
    """Apply the substitution. Returns the number of replacements made."""
    if not loc.path.exists():
        raise SystemExit(f"missing: {loc.path}")
    text = loc.path.read_text(encoding="utf-8")
    replacement = loc.template.format(version=new_version)
    new_text, n = loc.pattern.subn(replacement, text)
    if n == 0:
        raise SystemExit(
            f"version pattern not found in {loc.path} ({loc.label}). "
            f"Update bump_version.py to match the new layout."
        )
    if write and new_text != text:
        loc.path.write_text(new_text, encoding="utf-8")
    return n


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Bump the qPCR Analyzer version across the codebase.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "version",
        nargs="?",
        help="Explicit new version (e.g. 2.0.0). Mutually exclusive with --bump.",
    )
    ap.add_argument(
        "--bump",
        choices=["major", "minor", "patch"],
        help="Semver bump from current. Mutually exclusive with explicit version.",
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="Dry-run: show what would change, don't write files.",
    )
    args = ap.parse_args()

    current = read_current_version()

    if args.version is None and args.bump is None:
        print(f"current version: {current}")
        for loc in LOCATIONS:
            print(f"  - {loc.label}  [{loc.path.relative_to(ROOT)}]")
        return 0

    if args.version is not None and args.bump is not None:
        ap.error("pass either an explicit version or --bump, not both")

    if args.version is not None:
        if not SEMVER_FULL.match(args.version):
            ap.error(f"version must be X.Y.Z, got {args.version!r}")
        new_version = args.version
    else:
        new_version = bump(current, args.bump)

    if new_version == current:
        print(f"already at {current}; nothing to do")
        return 0

    print(f"bumping {current} -> {new_version}{'  (dry-run)' if args.check else ''}")
    for loc in LOCATIONS:
        n = update_location(loc, new_version, write=not args.check)
        print(f"  [{n}x] {loc.label}")
    if args.check:
        print("dry-run: no files modified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
