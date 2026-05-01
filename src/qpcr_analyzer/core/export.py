"""Excel and CSV-zip export of analysis results.

:func:`results_to_xlsx_bytes` returns an in-memory ``.xlsx`` workbook with one
``dCt_{HK}`` / ``ddCt_{HK}`` sheet per housekeeping gene plus matching
``formatted_*`` sheets — a wide grouped table (one block per target gene, one
column per biological group) with sample names included so the row provenance
is preserved. :func:`results_to_csv_zip_bytes` returns the same logical sheets
as separate CSV files inside a single zip archive.

The formatted-sheet builder runs a parity check at build time: every value
written to the formatted sheet is verified against the long-format result
DataFrame for the matching (target, group, sample) and a ``ValueError`` is
raised if any cell disagrees.
"""

from __future__ import annotations

import io
import re
import zipfile
from collections.abc import Iterator

import numpy as np
import pandas as pd

from .summary import sort_wells, target_order


def _safe_sheet_name(prefix: str, name: str) -> str:
    """Strip Excel-illegal characters and clamp to the 31-char sheet limit."""
    cleaned = re.sub(r"[\[\]*?\\/:']", "_", name)
    return (prefix + cleaned)[:31]


def _safe_filename(prefix: str, name: str) -> str:
    """Filesystem-safe filename for CSV-zip entries."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return f"{prefix}{cleaned}"


def _build_formatted_table(
    result_df: pd.DataFrame,
    value_col: str,
    value_label: str = "Relative Expression",
) -> pd.DataFrame:
    """Build a wide grouped table with explicit sample names per row.

    Layout: one block per target gene, separated by a blank row. Within each
    block, every biological group contributes two adjacent columns: one
    holding the sample name and one holding the value. Different groups can
    have different sample counts — short groups are padded with blanks.

    Column labels follow the pattern
    ``Target | Sample ({group1}) | {value_label} ({group1}) | ...``.

    A parity check verifies every emitted value matches the input
    ``result_df`` for the matching (Target, Group, Sample) cell, raising
    :class:`ValueError` on mismatch.
    """
    targets = list(result_df["Target"].unique())
    groups = list(result_df["Group"].unique())

    sample_col_for = {g: f"Sample ({g})" for g in groups}
    value_col_for = {g: f"{value_label} ({g})" for g in groups}
    columns = ["Target"]
    for g in groups:
        columns.extend([sample_col_for[g], value_col_for[g]])
    blank_row = {c: "" for c in columns}

    all_rows: list[dict] = []
    written: list[tuple[str, str, str, float | None]] = []  # for parity check
    for target in targets:
        sub = result_df[result_df["Target"] == target]
        header_row = {c: "" for c in columns}
        header_row["Target"] = target
        all_rows.append(header_row)

        per_group: dict[str, list[tuple[str, float | None]]] = {g: [] for g in groups}
        for g in groups:
            sub_g = sub.loc[sub["Group"] == g, ["Sample", value_col]]
            per_group[g] = [
                (str(r.Sample), None if pd.isna(getattr(r, value_col)) else float(getattr(r, value_col)))
                for r in sub_g.itertuples()
            ]
        max_len = max((len(v) for v in per_group.values()), default=0)
        for i in range(max_len):
            row: dict = {"Target": ""}
            for g in groups:
                if i < len(per_group[g]):
                    sample, value = per_group[g][i]
                    row[sample_col_for[g]] = sample
                    row[value_col_for[g]] = value
                    written.append((str(target), str(g), sample, value))
                else:
                    row[sample_col_for[g]] = ""
                    row[value_col_for[g]] = None
            all_rows.append(row)

        all_rows.append(dict(blank_row))

    out = pd.DataFrame(all_rows, columns=columns)
    _verify_formatted_against_source(written, result_df, value_col)
    return out


def _verify_formatted_against_source(
    written: list[tuple[str, str, str, float | None]],
    result_df: pd.DataFrame,
    value_col: str,
) -> None:
    """Cross-check every value placed in the formatted sheet vs. the source."""
    if not written:
        return
    source = result_df[["Target", "Group", "Sample", value_col]].copy()
    source["Target"] = source["Target"].astype(str)
    source["Group"] = source["Group"].astype(str)
    source["Sample"] = source["Sample"].astype(str)
    grouped = source.groupby(["Target", "Group", "Sample"], sort=False)[value_col].apply(list)

    used_indices: dict[tuple[str, str, str], int] = {}
    for target, group, sample, value in written:
        key = (target, group, sample)
        try:
            candidates = grouped.loc[key]
        except KeyError as ex:
            raise ValueError(
                f"Formatted-sheet parity check failed: "
                f"({target}, {group}, {sample}) not found in source data."
            ) from ex
        idx = used_indices.get(key, 0)
        used_indices[key] = idx + 1
        if idx >= len(candidates):
            raise ValueError(
                f"Formatted-sheet parity check failed: "
                f"({target}, {group}, {sample}) over-emitted "
                f"({idx + 1} > {len(candidates)})."
            )
        expected = candidates[idx]
        if pd.isna(expected) and value is None:
            continue
        if pd.isna(expected) or value is None:
            raise ValueError(
                f"Formatted-sheet parity check failed for "
                f"({target}, {group}, {sample}): NaN mismatch."
            )
        if not np.isclose(float(expected), float(value), atol=1e-9, rtol=1e-9):
            raise ValueError(
                f"Formatted-sheet parity check failed for "
                f"({target}, {group}, {sample}): "
                f"formatted={value} vs source={expected}."
            )


def _iter_sheets(
    raw_df: pd.DataFrame,
    dct_results: dict[str, pd.DataFrame],
    ddct_results: dict[str, pd.DataFrame],
) -> Iterator[tuple[str, str, pd.DataFrame]]:
    """Yield ``(prefix, name, dataframe)`` triples for every output sheet.

    Same logical content as the Excel workbook, used by both the xlsx and
    csv-zip writers so the two formats stay in sync. Prefixes are short
    keys for sheet/filename construction.
    """
    targets = target_order(raw_df) if "Target" in raw_df.columns else None
    if targets is not None:
        sorted_raw = sort_wells(raw_df, targets=targets)
    else:
        sorted_raw = raw_df
    yield ("raw_data", "", sorted_raw)

    for ref, df in dct_results.items():
        yield ("dCt_", ref, df)
        # Formatted ΔCt sheet displays relative expression vs HK
        # (2^−ΔCt), not the raw ΔCt values.
        yield (
            "formatted_dCt_",
            ref,
            _build_formatted_table(
                df, value_col="Expr_vs_HK", value_label="Relative Expression"
            ),
        )

    for ref, df in ddct_results.items():
        yield ("ddCt_", ref, df)
        yield (
            "formatted_ddCt_",
            ref,
            _build_formatted_table(
                df, value_col="Relative_Expr", value_label="Relative Expression"
            ),
        )


def results_to_xlsx_bytes(
    raw_df: pd.DataFrame,
    dct_results: dict[str, pd.DataFrame],
    ddct_results: dict[str, pd.DataFrame],
) -> bytes:
    """Write all results to an Excel workbook and return the bytes.

    Sheet layout
    ------------
    raw_data              standardised well-level data (with Excluded flag)
    dCt_{ref}             ΔCt table per housekeeping gene
    ddCt_{ref}            ΔΔCt table per housekeeping gene (includes Batch and
                          Reference_Group columns so provenance is self-documenting)
    formatted_dCt_{ref}   Wide grouped table with sample names per row
    formatted_ddCt_{ref}  Wide grouped table with sample names per row
    """
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        for prefix, name, df in _iter_sheets(raw_df, dct_results, ddct_results):
            sheet_name = (
                "raw_data" if prefix == "raw_data" else _safe_sheet_name(prefix, name)
            )
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    return buf.getvalue()


def results_to_csv_zip_bytes(
    raw_df: pd.DataFrame,
    dct_results: dict[str, pd.DataFrame],
    ddct_results: dict[str, pd.DataFrame],
) -> bytes:
    """Pack the same logical sheets as CSVs inside a single zip archive.

    Each Excel sheet maps to one ``.csv`` file inside the zip, named the same
    way as the Excel sheet (without the 31-char clamp).
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for prefix, name, df in _iter_sheets(raw_df, dct_results, ddct_results):
            filename = (
                "raw_data.csv"
                if prefix == "raw_data"
                else _safe_filename(prefix, name) + ".csv"
            )
            zf.writestr(filename, df.to_csv(index=False))
    return buf.getvalue()
