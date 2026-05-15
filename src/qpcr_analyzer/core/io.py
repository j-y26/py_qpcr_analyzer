"""File loaders for qPCR result tables.

Supports the file types instruments commonly export. Format dispatch is purely
extension-based — the ``filename`` argument is required even when ``source`` is
an in-memory stream so the extension can still be inspected.

Many qPCR exporters prepend a metadata block to the results table (Bio-Rad
CFX, Applied Biosystems QuantStudio, etc.) using a variety of prefixes —
``#``, ``*``, ``;``, ``//`` — or no prefix at all (e.g. ``Run Name: Plate1``).
The loader auto-detects the actual header row by scanning the first ~100
rows/lines and picking the earliest one whose cells contain qPCR-typical
keywords (well, target, sample, cq/ct, etc.). The same detection is applied
to Excel and CSV/TSV/TXT so the downstream column-mapping step sees a clean,
correctly-framed table.
"""

from __future__ import annotations

import csv
import io as _io
import re
from pathlib import Path
from typing import IO, Sequence, Union

import pandas as pd

SUPPORTED = {".xlsx", ".xls", ".csv", ".tsv", ".txt"}

# Lines starting with any of these are treated as metadata/comments and
# skipped when locating the header row.
_COMMENT_PREFIXES = ("#", "*", ";", "//", "!")

# Substrings (after lower-casing and stripping non-alphanumerics) that suggest
# a cell is part of a qPCR results-table header. Used to score candidate
# header rows.
_HEADER_KEYWORDS = (
    "well",
    "target",
    "sample",
    "cq",
    "ct",
    "gene",
    "assay",
    "position",
    "fluor",
    "detector",
    "group",
    "treatment",
    "condition",
    "batch",
    "run",
    "plate",
    "content",
    "biological",
)

_MIN_HEADER_HITS = 2
_MAX_SCAN_ROWS = 100


def read_table(source: Union[str, Path, IO[bytes]], filename: str) -> pd.DataFrame:
    """Read a qPCR results file into a :class:`pandas.DataFrame`.

    Args:
        source: A filesystem path or any binary file-like object (e.g.
            ``io.BytesIO`` from a NiceGUI upload). Pandas handles either.
        filename: Name (or path) used **only** to determine the file
            extension; needed when ``source`` is a stream that has no name.

    Returns:
        A :class:`pandas.DataFrame` with one row per well. The loader strips
        leading metadata/comment rows and auto-detects the real header but
        does **not** rename or coerce columns — see
        :func:`qpcr_analyzer.core.columns.apply_mapping`.

    Raises:
        ValueError: extension not in :data:`SUPPORTED`.
    """
    ext = Path(filename).suffix.lower()
    if ext not in SUPPORTED:
        raise ValueError(
            f"Unsupported file type: {ext}. Supported: {sorted(SUPPORTED)}"
        )
    if ext in {".xlsx", ".xls"}:
        return _read_excel(source)
    return _read_delimited(source, ext)


# ── Excel ────────────────────────────────────────────────────────────────────


def _read_excel(source: Union[str, Path, IO[bytes]]) -> pd.DataFrame:
    """Read .xlsx/.xls with header-row auto-detection.

    Reads the first sheet with ``header=None`` so prepended metadata rows
    are visible, then promotes the detected header row and drops everything
    above it. ``.xls`` requires the optional ``xlrd`` extra.
    """
    raw = pd.read_excel(source, header=None, dtype=object)
    rows = [
        [_cell_str(v) for v in raw.iloc[i].tolist()]
        for i in range(min(len(raw), _MAX_SCAN_ROWS))
    ]
    header_idx = _find_header_row(rows)
    headers = _dedup_headers([_cell_str(v) for v in raw.iloc[header_idx].tolist()])
    body = raw.iloc[header_idx + 1 :].copy()
    body.columns = headers
    return _post_clean(body)


# ── Delimited (csv / tsv / txt) ──────────────────────────────────────────────


def _read_delimited(
    source: Union[str, Path, IO[bytes]], ext: str
) -> pd.DataFrame:
    """Read CSV/TSV/TXT with header-row auto-detection and separator sniffing.

    Strategy: decode the bytes once, locate the header line by scanning for
    qPCR-keyword cells, then hand the trimmed text to ``pd.read_csv``. This
    handles metadata blocks regardless of whether they use ``#``, ``*``,
    ``;``, ``//`` or no prefix at all (e.g. ``Run Name: Plate1``).
    """
    raw_bytes = _read_bytes(source)
    # utf-8-sig strips a BOM if present (Excel-exported CSVs often carry one).
    text = raw_bytes.decode("utf-8-sig", errors="replace")
    lines = text.splitlines()

    if ext == ".csv":
        sep = ","
    elif ext == ".tsv":
        sep = "\t"
    else:  # .txt
        sep = _sniff_separator(lines)

    # Parse just enough rows to locate the header, using the stdlib CSV
    # reader so quoting and embedded separators are honoured.
    sample = "\n".join(lines[:_MAX_SCAN_ROWS])
    try:
        scan_rows = list(csv.reader(_io.StringIO(sample), delimiter=sep))
    except csv.Error:
        scan_rows = [ln.split(sep) for ln in lines[:_MAX_SCAN_ROWS]]
    header_idx = _find_header_row(scan_rows)

    body_text = "\n".join(lines[header_idx:])
    df = pd.read_csv(
        _io.StringIO(body_text),
        sep=sep,
        comment="#",
        skip_blank_lines=True,
        engine="python",
    )
    return _post_clean(df)


def _sniff_separator(lines: Sequence[str]) -> str:
    """Guess the column separator for a .txt file from its first few lines.

    Picks whichever of ``\\t`` / ``,`` / ``;`` gives the most consistent
    non-trivial field count across non-comment lines. Falls back to TSV
    because qPCR ``.txt`` exports are tab-separated more often than not.
    """
    candidates = ["\t", ",", ";"]
    best = ("\t", 0)
    for sep in candidates:
        counts = []
        for ln in lines[:30]:
            stripped = ln.strip()
            if not stripped or any(
                stripped.startswith(p) for p in _COMMENT_PREFIXES
            ):
                continue
            n = len(ln.split(sep))
            if n > 1:
                counts.append(n)
        if not counts:
            continue
        # Score = most-common count × number of lines that hit it.
        mode_count = max(set(counts), key=counts.count)
        score = mode_count * counts.count(mode_count)
        if score > best[1]:
            best = (sep, score)
    return best[0]


# ── Header detection ─────────────────────────────────────────────────────────


def _norm_token(s: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _cell_str(v: object) -> str:
    """Render a raw cell value to a stripped string, mapping NaN to ``''``."""
    if v is None:
        return ""
    # pandas missing values
    try:
        if pd.isna(v):  # handles NaN, NaT, pd.NA
            return ""
    except (TypeError, ValueError):
        pass
    return str(v).strip()


def _row_keyword_score(fields: Sequence[object]) -> int:
    """Number of cells whose normalised text contains a known header keyword."""
    score = 0
    for f in fields:
        n = _norm_token(f)
        if not n:
            continue
        if any(k in n for k in _HEADER_KEYWORDS):
            score += 1
    return score


def _find_header_row(rows: Sequence[Sequence[object]]) -> int:
    """Return the index of the row most likely to be the column header.

    Picks the earliest row (within ``rows``) with the highest count of cells
    containing qPCR-typical keywords (well/target/sample/cq/ct/…). Rows
    whose first cell starts with a comment prefix (``#``, ``*``, ``;``,
    ``//``, ``!``) and rows that are entirely blank are skipped. Falls back
    to row 0 if no row scores at least :data:`_MIN_HEADER_HITS`.
    """
    best_score = 0
    best_idx = 0
    for i, fields in enumerate(rows):
        cells = [_cell_str(f) for f in fields]
        if not any(cells):
            continue
        first = next((c for c in cells if c), "")
        if any(first.startswith(p) for p in _COMMENT_PREFIXES):
            continue
        score = _row_keyword_score(cells)
        if score > best_score:
            best_score = score
            best_idx = i
    return best_idx if best_score >= _MIN_HEADER_HITS else 0


def _dedup_headers(headers: Sequence[str]) -> list[str]:
    """Replace blanks with ``Column_N`` and disambiguate duplicates."""
    seen: dict[str, int] = {}
    out: list[str] = []
    for i, h in enumerate(headers):
        name = h.strip() if h else ""
        if not name or name.lower() == "nan":
            name = f"Column_{i + 1}"
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        out.append(name)
    return out


# ── Post-load cleanup ────────────────────────────────────────────────────────


def _post_clean(df: pd.DataFrame) -> pd.DataFrame:
    """Drop fully-empty rows and trim columns whose header is blank.

    Row-level cleanup of missing Target/Sample happens later in
    :func:`qpcr_analyzer.core.columns.apply_mapping` once roles are known.
    """
    df = df.dropna(axis=0, how="all").reset_index(drop=True)
    # Drop columns whose header is blank/NaN AND whose values are all blank.
    keep_cols = []
    for c in df.columns:
        name = str(c).strip()
        if name and not name.lower().startswith("unnamed"):
            keep_cols.append(c)
            continue
        if df[c].notna().any() and df[c].astype(str).str.strip().ne("").any():
            keep_cols.append(c)
    return df[keep_cols]


# ── Helpers ──────────────────────────────────────────────────────────────────


def _read_bytes(source: Union[str, Path, IO[bytes]]) -> bytes:
    """Slurp the input source into bytes, leaving the caller's stream intact."""
    if isinstance(source, (str, Path)):
        return Path(source).read_bytes()
    pos = source.tell() if hasattr(source, "tell") else None
    data = source.read()
    if pos is not None and hasattr(source, "seek"):
        source.seek(pos)
    if isinstance(data, str):
        return data.encode("utf-8")
    return data
