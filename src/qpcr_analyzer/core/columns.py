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

ROLES = ("well", "target", "sample", "cq", "group")
ROLE_LABELS = {
    "well": "Well",
    "target": "Target",
    "sample": "Sample",
    "cq": "Cq",
    "group": "Group",
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


def apply_mapping(df: pd.DataFrame, mapping: ColumnMapping) -> pd.DataFrame:
    """Project ``df`` to canonical columns and coerce types.

    Returns a new DataFrame containing only the mapped columns, renamed to
    ``Well``, ``Target``, ``Sample``, ``Cq`` and ``Group``. Behaviour:

    * If the ``group`` role is unmapped, every row is assigned to ``"all"``
      so downstream code can rely on ``Group`` always being present.
    * Missing/null values in ``Group`` become ``"unassigned"``.
    * ``Cq`` is coerced to numeric (non-numeric values become ``NaN`` and
      will be flagged as outliers downstream).
    * ``Well`` / ``Target`` / ``Sample`` are forced to string.
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
    out["Group"] = out["Group"].fillna("unassigned").astype(str)
    out["Cq"] = pd.to_numeric(out["Cq"], errors="coerce")
    out["Well"] = out["Well"].astype(str)
    out["Target"] = out["Target"].astype(str)
    out["Sample"] = out["Sample"].astype(str)
    return out
