"""Excel export of analysis results.

:func:`results_to_xlsx_bytes` produces an in-memory ``.xlsx`` workbook the
NiceGUI app streams to the browser as a download. The workbook intentionally
mirrors the in-memory result dicts: one ``dCt_{HK}`` and ``ddCt_{HK}`` sheet
per housekeeping gene, plus matching ``prism_*`` sheets pre-formatted for
direct copy-paste into a GraphPad Prism *Grouped* table.
"""

from __future__ import annotations

import io
import re

import pandas as pd


def _safe_sheet_name(prefix: str, name: str) -> str:
    """Strip Excel-illegal characters and clamp to the 31-char sheet limit."""
    cleaned = re.sub(r"[\[\]*?\\/:']", "_", name)
    return (prefix + cleaned)[:31]


def _build_prism_table(
    result_df: pd.DataFrame,
    value_col: str,
) -> pd.DataFrame:
    """Build a GraphPad Prism-ready grouped table.

    Layout: one block per target gene, separated by a blank row.  Within each
    block the columns are biological groups and the rows are per-sample values
    (replicate stacking).  Suitable for copy-paste into a Prism Grouped table.

    Example (Relative_Expr, two targets, two groups)::

        Target         Control    Treatment
        GeneX
                       1.05       2.31
                       0.93       2.15
                       1.02       2.48

        GeneY
                       0.88       1.72
                       ...
    """
    targets = list(result_df["Target"].unique())
    groups = list(result_df["Group"].unique())

    all_rows: list[dict] = []
    for target in targets:
        sub = result_df[result_df["Target"] == target]
        # Header row: target name in first col, group names in value cols
        all_rows.append({"Target": target, **{g: g for g in groups}})
        # Data rows: one per sample
        per_group: dict[str, list] = {g: [] for g in groups}
        for g in groups:
            per_group[g] = sub.loc[sub["Group"] == g, value_col].tolist()
        max_len = max((len(v) for v in per_group.values()), default=0)
        for i in range(max_len):
            all_rows.append({
                "Target": "",
                **{g: (per_group[g][i] if i < len(per_group[g]) else None) for g in groups},
            })
        # Blank separator
        all_rows.append({"Target": "", **{g: None for g in groups}})

    return pd.DataFrame(all_rows)


def results_to_xlsx_bytes(
    raw_df: pd.DataFrame,
    dct_results: dict[str, pd.DataFrame],
    ddct_results: dict[str, pd.DataFrame],
) -> bytes:
    """Write all results to an Excel workbook and return the bytes.

    Sheet layout
    ------------
    raw_data          standardised well-level data (with Excluded flag)
    dCt_{ref}         ΔCt table per housekeeping gene
    ddCt_{ref}        ΔΔCt table per housekeeping gene (includes Batch and
                      Reference_Group columns so provenance is self-documenting)
    prism_dCt_{ref}   Prism-ready grouped table of dCt values
    prism_ddCt_{ref}  Prism-ready grouped table of Relative_Expr values
    """
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        raw_df.to_excel(writer, sheet_name="raw_data", index=False)

        for ref, df in dct_results.items():
            df.to_excel(
                writer, sheet_name=_safe_sheet_name("dCt_", ref), index=False
            )
            prism = _build_prism_table(df, value_col="dCt")
            prism.to_excel(
                writer, sheet_name=_safe_sheet_name("prism_dCt_", ref), index=False
            )

        for ref, df in ddct_results.items():
            df.to_excel(
                writer, sheet_name=_safe_sheet_name("ddCt_", ref), index=False
            )
            prism = _build_prism_table(df, value_col="Relative_Expr")
            prism.to_excel(
                writer, sheet_name=_safe_sheet_name("prism_ddCt_", ref), index=False
            )

    return buf.getvalue()
