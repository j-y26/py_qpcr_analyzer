"""Column-role detection, mapping, and sample/group validation.

Different qPCR instruments use different column names (Applied Biosystems'
``Well Position`` / ``Target Name`` / ``CT``, Bio-Rad CFX's ``Well`` /
``Target`` / ``Cq``, etc.). This module reduces them to a canonical set of
five roles — ``well``, ``target``, ``sample``, ``cq``, ``group`` — using
exact, substring, and fuzzy matching, with two heuristic fallbacks for files
whose headers don't match any synonym.

Public API:
    :data:`ROLES`, :data:`ROLE_LABELS`
        Ordered role tuple and their display labels.
    :class:`ColumnMapping`
        Result object with ``assignments`` and per-role confidence.
    :func:`detect_columns`
        Score every column for every role, pick the best, run heuristics.
    :func:`apply_mapping`
        Rename detected columns to canonical names and coerce types.
    :func:`validate_sample_groups`
        Cross-check that each sample belongs to exactly one group.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd
from rapidfuzz import fuzz

ROLES = ("well", "target", "sample", "cq", "group", "batch")
ROLE_LABELS = {
    "well": "Well",
    "target": "Target",
    "sample": "Sample",
    "cq": "Cq",
    "group": "Group",
    "batch": "Batch",
}
REQUIRED = ("well", "target", "sample", "cq")

SYNONYMS: dict[str, list[str]] = {
    "well": ["well", "wells", "wellid", "location", "wellposition", "position"],
    "target": ["target", "gene", "assay", "assayname", "targetname", "detector"],
    "sample": ["sample", "sampleid", "samplename"],
    "cq": ["cq", "ct", "cqvalue", "ctvalue", "cqmean", "ctmean"],
    "group": [
        "group",
        "groupid",
        "groupname",
        "condition",
        "treatment",
        "biologicalsetname",
    ],
    "batch": [
        "batch",
        "batches",
        "batchid",
        "batchname",
        "run",
        "runid",
        "runname",
        "plate",
        "plateid",
        "platename",
        "experiment",
        "experimentid",
    ],
}

WELL_RE = re.compile(r"^[A-Pa-p]\d{1,2}$")
ACCEPT_THRESHOLD = 0.85


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _has_well_pattern(series: pd.Series, min_frac: float = 0.8) -> bool:
    """Return True if >= min_frac of non-null values match A1-style well IDs."""
    vals = series.dropna().astype(str)
    if len(vals) == 0:
        return False
    return float(vals.str.match(WELL_RE).mean()) >= min_frac


@dataclass
class ColumnMapping:
    """Result of :func:`detect_columns`.

    Attributes:
        assignments: ``{role: original_column_name | None}``. Roles missing
            from the input have value ``None``.
        confidence: ``{role: float in [0, 1]}``, the score of the picked
            column. Used by the UI to display "auto · 95 %" badges.
    """

    assignments: dict[str, str | None] = field(default_factory=dict)
    confidence: dict[str, float] = field(default_factory=dict)

    def validate(self) -> list[str]:
        """Return a list of error strings; empty when the mapping is usable.

        Checks that all four required roles (``well``, ``target``,
        ``sample``, ``cq``) are filled and that no two roles point at the
        same source column.
        """
        errors: list[str] = []
        for role in REQUIRED:
            if not self.assignments.get(role):
                errors.append(f"Missing required column for {ROLE_LABELS[role]}")
        assigned = [v for v in self.assignments.values() if v]
        if len(assigned) != len(set(assigned)):
            errors.append("The same column is mapped to multiple roles")
        return errors


def detect_columns(df: pd.DataFrame) -> ColumnMapping:
    """Detect which DataFrame column plays each role.

    Algorithm:
        1. Normalise every header (lowercase, strip non-alphanumerics).
        2. For every (column, role) pair compute a score:
           exact-match → 1.0, substring → 0.95, contained → 0.9,
           else fuzzy ``rapidfuzz.fuzz.ratio / 100``.
        3. For each role, pick the highest-scoring column above
           :data:`ACCEPT_THRESHOLD`. On a tie for the *well* role, prefer
           a column whose values look like A1-format wells — handles
           Applied Biosystems output that has both ``Well`` (integers) and
           ``Well Position`` (A1, A2, …).
        4. Heuristic fallbacks for unassigned columns: a numeric column
           with median in [5, 45] is assumed to be Cq; an alphanumeric
           A1-pattern column is assumed to be Well.
    """
    cols = list(df.columns)
    norm_cols = {c: _norm(c) for c in cols}
    mapping = ColumnMapping()

    # Score every column for every role; track all candidates at the top score.
    all_scored: dict[str, list[tuple[float, str]]] = {}
    for role in ROLES:
        syns = SYNONYMS[role]
        scored: list[tuple[float, str]] = []
        for c, nc in norm_cols.items():
            if nc in syns:
                score = 1.0
            elif any(s in nc for s in syns):
                score = 0.95
            elif any(nc in s and len(nc) >= 2 for s in syns):
                score = 0.9
            else:
                score = max((fuzz.ratio(nc, s) for s in syns), default=0) / 100.0
            scored.append((score, c))
        scored.sort(key=lambda x: -x[0])
        all_scored[role] = scored

    for role in ROLES:
        scored = all_scored[role]
        if not scored or scored[0][0] < ACCEPT_THRESHOLD:
            mapping.assignments[role] = None
            mapping.confidence[role] = scored[0][0] if scored else 0.0
            continue

        best_score = scored[0][0]
        # All candidates tied at the best score
        tied = [c for s, c in scored if s >= best_score - 1e-9]

        if role == "well" and len(tied) > 1:
            # When multiple columns score equally (e.g. "Well" and "Well Position"
            # on Applied Biosystems output), prefer the one whose values actually
            # contain A1-format well IDs rather than bare integers.
            well_cols = [c for c in tied if _has_well_pattern(df[c])]
            best_col = well_cols[0] if well_cols else tied[0]
        else:
            best_col = tied[0]

        mapping.assignments[role] = best_col
        mapping.confidence[role] = best_score

    assigned = {v for v in mapping.assignments.values() if v}
    unassigned_cols = [c for c in cols if c not in assigned]

    # Heuristic fallback: numeric column in qPCR Cq range → Cq
    if not mapping.assignments.get("cq") and unassigned_cols:
        for c in list(unassigned_cols):
            numeric = pd.to_numeric(df[c], errors="coerce")
            valid = numeric.dropna()
            if len(valid) >= max(3, int(len(df) * 0.5)) and 5 <= valid.median() <= 45:
                mapping.assignments["cq"] = c
                mapping.confidence["cq"] = max(mapping.confidence.get("cq", 0), 0.6)
                assigned.add(c)
                unassigned_cols.remove(c)
                break

    # Heuristic fallback: column whose values look like well IDs → Well
    if not mapping.assignments.get("well") and unassigned_cols:
        for c in list(unassigned_cols):
            if _has_well_pattern(df[c]):
                mapping.assignments["well"] = c
                mapping.confidence["well"] = max(mapping.confidence.get("well", 0), 0.6)
                assigned.add(c)
                unassigned_cols.remove(c)
                break

    return mapping


def validate_sample_groups(df: pd.DataFrame) -> list[str]:
    """Check that each sample belongs to exactly one group.

    Args:
        df: DataFrame with 'Sample' and 'Group' columns (standardised names).

    Returns:
        List of error strings, empty when data is valid.
    """
    errors: list[str] = []
    per_sample = df.groupby("Sample")["Group"].nunique()
    multi = per_sample[per_sample > 1].index.tolist()
    for sample in multi:
        groups = sorted(df.loc[df["Sample"] == sample, "Group"].unique())
        errors.append(f"Sample '{sample}' belongs to multiple groups: {groups}")
    return errors


def validate_sample_batches(df: pd.DataFrame) -> list[str]:
    """Check that each sample belongs to exactly one batch.

    Args:
        df: DataFrame with 'Sample' and 'Batch' columns (standardised names).
            If 'Batch' is absent, returns an empty list (nothing to validate).

    Returns:
        List of error strings, empty when data is valid.
    """
    if "Batch" not in df.columns:
        return []
    errors: list[str] = []
    per_sample = df.groupby("Sample")["Batch"].nunique()
    multi = per_sample[per_sample > 1].index.tolist()
    for sample in multi:
        batches = sorted(df.loc[df["Sample"] == sample, "Batch"].unique())
        errors.append(f"Sample '{sample}' belongs to multiple batches: {batches}")
    return errors


# Strings users tend to type for "no signal" in a qPCR exporter — coerced to
# NaN before the numeric conversion so well-typed numbers in the same column
# aren't dragged to NaN by ``pd.to_numeric``'s all-or-nothing dtype inference.
_CQ_MISSING_TOKENS = frozenset(
    {
        "",
        "na",
        "n/a",
        "nan",
        "null",
        "none",
        "nd",
        "n.d.",
        "-",
        "--",
        "undetermined",
        "undet",
        "no ct",
        "no cq",
    }
)


def _coerce_cq(series: pd.Series) -> pd.Series:
    """Coerce a Cq column to numeric, tolerating whitespace and "no signal" tokens.

    ``pd.to_numeric(errors='coerce')`` already maps unparseable strings to
    NaN, but qPCR exporters tend to write tokens like ``Undetermined`` /
    ``N/A`` interleaved with numeric Cq values. We strip whitespace and map
    those tokens explicitly so a stray trailing space or unusual exporter
    wording cannot poison the dtype inference.
    """
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    s = series.astype(str).str.strip()
    lowered = s.str.lower()
    s = s.mask(lowered.isin(_CQ_MISSING_TOKENS))
    return pd.to_numeric(s, errors="coerce")


def _coerce_label(series: pd.Series) -> pd.Series:
    """Coerce a label column (Well/Target/Sample) to stripped string.

    NaN-like values become empty strings so callers can test ``== ""`` to
    spot incomplete annotations.
    """
    s = series.astype("object").where(series.notna(), "")
    s = s.astype(str).str.strip()
    return s.replace({"nan": "", "None": "", "NaN": ""})


def apply_mapping(df: pd.DataFrame, mapping: ColumnMapping) -> pd.DataFrame:
    """Project ``df`` to canonical columns, coerce types, drop unusable rows.

    Returns a new DataFrame containing only the mapped columns, renamed to
    ``Well``, ``Target``, ``Sample``, ``Cq`` and ``Group``. Behaviour:

    * If the ``group`` role is unmapped, every row is assigned to ``"all"``
      so downstream code can rely on ``Group`` always being present.
    * Missing/null values in ``Group`` become ``"unassigned"``.
    * ``Cq`` is coerced to numeric (whitespace, ``Undetermined``, ``N/A``
      and similar markers become ``NaN`` and will be flagged as outliers
      downstream).
    * ``Well`` / ``Target`` / ``Sample`` are forced to stripped string.
    * Rows where **all of** Cq, Sample, and Target are blank are dropped
      silently — these are truly unused wells the exporter wrote out with
      no data. Rows missing only *some* required fields are **kept** so
      :func:`find_incomplete_rows` (and the UI) can surface them for user
      review rather than silently discarding partial annotations.
    """
    rename: dict[str, str] = {}
    for role, col in mapping.assignments.items():
        if col is None:
            continue
        rename[col] = ROLE_LABELS[role]
    keep = list(rename.keys())
    out = df[keep].rename(columns=rename).copy()
    if "Group" not in out.columns:
        out["Group"] = "all"
    out["Group"] = out["Group"].fillna("unassigned").astype(str).str.strip()
    out.loc[out["Group"] == "", "Group"] = "unassigned"
    out["Cq"] = _coerce_cq(out["Cq"])
    out["Well"] = _coerce_label(out["Well"])
    out["Target"] = _coerce_label(out["Target"])
    out["Sample"] = _coerce_label(out["Sample"])
    if "Batch" in out.columns:
        out["Batch"] = out["Batch"].fillna("batch_1").astype(str).str.strip()
        out.loc[out["Batch"] == "", "Batch"] = "batch_1"
    # Silent drop: only rows that have NO biology data at all. Rows with
    # partial missingness (e.g. Sample missing but Target+Cq present) stay
    # in the frame so the UI can flag them; the UI is responsible for
    # excluding them with explicit user acknowledgement.
    fully_empty = (
        out["Cq"].isna() & (out["Target"] == "") & (out["Sample"] == "")
    )
    return out.loc[~fully_empty].reset_index(drop=True)


# Fields a row must have annotated to be analyzable downstream. Cq is not
# in this list — a row with valid Target+Sample but NaN Cq is an
# Undetermined well, biologically meaningful and tracked through the
# outlier-review UI.
_ANNOTATION_FIELDS: tuple[str, ...] = ("Well", "Target", "Sample")


def find_incomplete_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows with one or more required annotations missing.

    "Incomplete" here means at least one of Well / Target / Sample is blank
    or NaN. Cq missing on its own is **not** flagged: it represents an
    Undetermined well, valid biology surfaced through outlier review.

    Args:
        df: Standardised DataFrame from :func:`apply_mapping`.

    Returns:
        A subset of ``df`` containing only the flagged rows, with an extra
        ``Missing`` column listing the names of the missing required fields
        (comma-joined, e.g. ``"Target"`` or ``"Sample, Well"``). The
        original row index is preserved so callers can map back to the
        source row order. Empty DataFrame when nothing is missing.
    """
    if df.empty:
        return df.assign(Missing=pd.Series(dtype=str)).iloc[0:0]

    masks: dict[str, pd.Series] = {}
    for field in _ANNOTATION_FIELDS:
        if field not in df.columns:
            continue
        col = df[field]
        if pd.api.types.is_numeric_dtype(col):
            masks[field] = col.isna()
        else:
            masks[field] = col.isna() | (col.astype(str).str.strip() == "")

    any_missing = pd.Series(False, index=df.index)
    for m in masks.values():
        any_missing = any_missing | m

    if not any_missing.any():
        return df.loc[any_missing].assign(Missing="")

    out = df.loc[any_missing].copy()
    def _missing_for_row(idx: int) -> str:
        return ", ".join(f for f, m in masks.items() if bool(m.loc[idx]))
    out["Missing"] = [_missing_for_row(i) for i in out.index]
    return out


def annotation_coverage(df: pd.DataFrame) -> dict[str, tuple[int, int]]:
    """Per required-annotation coverage: ``{field: (n_filled, n_total)}``.

    Used by the UI to render "Sample: 45 / 48 (94%)" style coverage badges
    in the column-mapping step. ``Cq`` is included alongside the annotation
    fields so the UI can show how many rows have a numeric Cq, but a low
    Cq coverage is informational only — it does not make a row "incomplete".
    """
    out: dict[str, tuple[int, int]] = {}
    n_total = len(df)
    for field in (*_ANNOTATION_FIELDS, "Cq"):
        if field not in df.columns:
            continue
        col = df[field]
        if pd.api.types.is_numeric_dtype(col):
            n_filled = int(col.notna().sum())
        else:
            n_filled = int(
                (col.notna() & (col.astype(str).str.strip() != "")).sum()
            )
        out[field] = (n_filled, n_total)
    return out


def drop_incomplete_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` with all rows flagged by :func:`find_incomplete_rows` removed.

    Convenience for the UI: after the user acknowledges the incomplete-rows
    table in step 2, this is what strips them out before the standardised
    DataFrame propagates to downstream steps.
    """
    incomplete = find_incomplete_rows(df)
    if incomplete.empty:
        return df.reset_index(drop=True)
    return df.drop(incomplete.index).reset_index(drop=True)
