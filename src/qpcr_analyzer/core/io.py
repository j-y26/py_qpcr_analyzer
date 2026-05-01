"""File loaders for qPCR result tables.

Supports the file types instruments commonly export. Format dispatch is purely
extension-based — the ``filename`` argument is required even when ``source`` is
an in-memory stream so the extension can still be inspected.
"""

from __future__ import annotations

from pathlib import Path
from typing import IO, Union

import pandas as pd

SUPPORTED = {".xlsx", ".xls", ".csv", ".tsv", ".txt"}


def read_table(source: Union[str, Path, IO[bytes]], filename: str) -> pd.DataFrame:
    """Read a qPCR results file into a :class:`pandas.DataFrame`.

    Args:
        source: A filesystem path or any binary file-like object (e.g.
            ``io.BytesIO`` from a NiceGUI upload). Pandas handles either.
        filename: Name (or path) used **only** to determine the file
            extension; needed when ``source`` is a stream that has no name.

    Returns:
        A :class:`pandas.DataFrame` with one row per well. No column renaming
        or type coercion is done here — see :func:`qpcr_analyzer.core.columns.apply_mapping`.

    Raises:
        ValueError: extension not in :data:`SUPPORTED`.
    """
    ext = Path(filename).suffix.lower()
    if ext not in SUPPORTED:
        raise ValueError(
            f"Unsupported file type: {ext}. Supported: {sorted(SUPPORTED)}"
        )
    if ext in {".xlsx", ".xls"}:
        # .xls requires the optional `xlrd` extra (`pip install qpcr-analyzer[xls]`).
        return pd.read_excel(source)
    # ``comment="#"`` skips the metadata header that Thermo Fisher Design &
    # Analysis prepends to its CSV/TSV exports (e.g. ``# File Name: ...``,
    # ``# Instrument Type: ...`` — about 25 lines before the real header).
    if ext == ".csv":
        return pd.read_csv(source, comment="#")
    return pd.read_csv(source, sep="\t", comment="#")

