from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

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
        "batch",
    ],
}

WELL_RE = re.compile(r"^[A-Pa-p]\d{1,2}$")
ACCEPT_THRESHOLD = 0.85


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


@dataclass
class ColumnMapping:
    assignments: dict[str, str | None] = field(default_factory=dict)
    confidence: dict[str, float] = field(default_factory=dict)

    def validate(self) -> list[str]:
        errors: list[str] = []
        for role in REQUIRED:
            if not self.assignments.get(role):
                errors.append(f"Missing required column for {ROLE_LABELS[role]}")
        assigned = [v for v in self.assignments.values() if v]
        if len(assigned) != len(set(assigned)):
            errors.append("The same column is mapped to multiple roles")
        return errors


def detect_columns(df: pd.DataFrame) -> ColumnMapping:
    cols = list(df.columns)
    norm_cols = {c: _norm(c) for c in cols}
    mapping = ColumnMapping()

    for role in ROLES:
        syns = SYNONYMS[role]
        best_col, best_score = None, 0.0
        for c, nc in norm_cols.items():
            if nc in syns:
                score = 1.0
            elif any(s in nc for s in syns):
                score = 0.95
            elif any(nc in s and len(nc) >= 2 for s in syns):
                score = 0.9
            else:
                score = max((fuzz.ratio(nc, s) for s in syns), default=0) / 100.0
            if score > best_score:
                best_col, best_score = c, score
        if best_score >= ACCEPT_THRESHOLD:
            mapping.assignments[role] = best_col
        else:
            mapping.assignments[role] = None
        mapping.confidence[role] = best_score

    assigned = {v for v in mapping.assignments.values() if v}
    unassigned_cols = [c for c in cols if c not in assigned]

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

    if not mapping.assignments.get("well") and unassigned_cols:
        for c in list(unassigned_cols):
            vals = df[c].dropna().astype(str)
            if len(vals) == 0:
                continue
            hits = vals.str.match(WELL_RE).sum()
            if hits / len(vals) > 0.8:
                mapping.assignments["well"] = c
                mapping.confidence["well"] = max(mapping.confidence.get("well", 0), 0.6)
                assigned.add(c)
                unassigned_cols.remove(c)
                break

    return mapping


def apply_mapping(df: pd.DataFrame, mapping: ColumnMapping) -> pd.DataFrame:
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
