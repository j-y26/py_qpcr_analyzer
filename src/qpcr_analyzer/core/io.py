from __future__ import annotations

from pathlib import Path
from typing import IO, Union

import pandas as pd

SUPPORTED = {".xlsx", ".xls", ".csv", ".tsv", ".txt"}


def read_table(source: Union[str, Path, IO[bytes]], filename: str) -> pd.DataFrame:
    ext = Path(filename).suffix.lower()
    if ext not in SUPPORTED:
        raise ValueError(
            f"Unsupported file type: {ext}. Supported: {sorted(SUPPORTED)}"
        )
    if ext in {".xlsx", ".xls"}:
        return pd.read_excel(source)
    if ext == ".csv":
        return pd.read_csv(source)
    return pd.read_csv(source, sep="\t")
