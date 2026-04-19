from __future__ import annotations

import io
import re

import pandas as pd


def _safe_sheet_name(name: str, prefix: str = "ref_") -> str:
    cleaned = re.sub(r"[\[\]\*\?\\/:]", "_", name)
    return (prefix + cleaned)[:31]


def results_to_xlsx_bytes(
    raw_df: pd.DataFrame,
    results: dict[str, pd.DataFrame],
) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        raw_df.to_excel(writer, sheet_name="raw_data", index=False)
        for ref, res in results.items():
            res.to_excel(writer, sheet_name=_safe_sheet_name(ref), index=False)
    return buf.getvalue()
