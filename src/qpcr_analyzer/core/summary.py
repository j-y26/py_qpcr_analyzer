"""Dataset summarisation and well-sort helpers.

Pure-Python utilities used by the UI to render the upload-time dataset
summary and to sort wells in a biology-friendly order: by *target* (in the
order targets first appear in the source file) → *group* → *sample* →
*well* (natural A1, A2 … A10 order, not lexicographic).

The helpers also bundle (Sample × Target) replicates into "blocks" for the
exclusion review panel, so the UI can show every replicate of a flagged
block together rather than just the individual flagged wells.
"""

from __future__ import annotations

import re
from datetime import datetime

import pandas as pd

WELL_RE = re.compile(r"^([A-Pa-p])(\d{1,2})$")


def well_sort_key(well: str) -> tuple:
    """Key for natural sorting of A1-format wells (A1, A2, ..., A10, A11)."""
    m = WELL_RE.match(str(well))
    if m:
        return (0, m.group(1).upper(), int(m.group(2)))
    return (1, str(well), 0)


def target_order(df: pd.DataFrame) -> list[str]:
    """Targets in the order they first appear in the source DataFrame."""
    return list(pd.unique(df["Target"].astype(str)))


def sample_order(df: pd.DataFrame) -> list[str]:
    """Sample names in the order they first appear in the source DataFrame."""
    return list(pd.unique(df["Sample"].astype(str)))


def group_order(df: pd.DataFrame) -> list[str]:
    """Group names in the order they first appear in the source DataFrame."""
    if "Group" not in df.columns:
        return []
    return list(pd.unique(df["Group"].astype(str)))


def sort_wells(
    df: pd.DataFrame,
    targets: list[str] | None = None,
) -> pd.DataFrame:
    """Sort by Target → Group → Sample → Well, all in *file-appearance* order.

    Replicate wells of the same (Sample × Target) end up adjacent. None of
    Group / Sample is sorted alphanumerically — biological labels often
    encode meaning (e.g. ``donor_3`` before ``donor_10``) that
    lexicographic sort would break, so we honour the order they first
    appear in the upload. Wells inside a block use natural A1, A2, …
    ordering since well IDs are just plate coordinates.

    Args:
        df: Standardised DataFrame with columns Well, Target, Sample, Group.
        targets: Explicit target order. When None, derived from ``df`` via
            :func:`target_order`.
    """
    if targets is None:
        targets = target_order(df)
    target_rank = {t: i for i, t in enumerate(targets)}
    group_rank = {g: i for i, g in enumerate(group_order(df))}
    sample_rank = {s: i for i, s in enumerate(sample_order(df))}

    out = df.copy()
    out["__t_rank"] = out["Target"].astype(str).map(target_rank).fillna(len(targets))
    out["__g_rank"] = (
        out["Group"].astype(str).map(group_rank).fillna(len(group_rank))
        if "Group" in out.columns
        else 0
    )
    out["__s_rank"] = (
        out["Sample"].astype(str).map(sample_rank).fillna(len(sample_rank))
    )
    out["__well_key"] = out["Well"].astype(str).map(well_sort_key)
    out = out.sort_values(
        by=["__t_rank", "__g_rank", "__s_rank", "__well_key"],
        kind="mergesort",
    ).drop(columns=["__t_rank", "__g_rank", "__s_rank", "__well_key"])
    return out.reset_index(drop=True)


def build_blocks(
    flagged: pd.DataFrame,
    excluded_wells: set[str],
    targets: list[str] | None = None,
) -> list[dict]:
    """Group wells into (Sample × Target) blocks, exclusion-bearing first.

    Each block is a dict with keys:

        target, group, sample, replicates (list of well-row dicts),
        n_excluded (int), n_flagged (int), has_exclusion (bool)

    Blocks where any replicate is in ``excluded_wells`` (or is auto-flagged
    NaN/outlier) come first, sorted by the same Target → Group → Sample key
    as :func:`sort_wells`. Clean blocks come after.
    """
    if targets is None:
        targets = target_order(flagged)
    target_rank = {t: i for i, t in enumerate(targets)}
    group_rank = {g: i for i, g in enumerate(group_order(flagged))}
    sample_rank = {s: i for i, s in enumerate(sample_order(flagged))}

    sorted_df = sort_wells(flagged, targets=targets)
    blocks: list[dict] = []
    grouped = sorted_df.groupby(["Target", "Group", "Sample"], sort=False)
    for (target, group, sample), sub in grouped:
        replicates = []
        n_excluded = 0
        n_flagged = 0
        for r in sub.itertuples():
            cq = None if pd.isna(r.Cq) else float(r.Cq)
            is_excluded = r.Well in excluded_wells
            is_flagged = bool(getattr(r, "Outlier", False))
            is_nan = cq is None
            if is_excluded:
                n_excluded += 1
            if is_flagged or is_nan:
                n_flagged += 1
            replicates.append(
                {
                    "well": str(r.Well),
                    "cq": cq,
                    "excluded": is_excluded,
                    "outlier": is_flagged,
                    "is_nan": is_nan,
                }
            )
        blocks.append(
            {
                "target": str(target),
                "group": str(group),
                "sample": str(sample),
                "replicates": replicates,
                "n_excluded": n_excluded,
                "n_flagged": n_flagged,
                "has_exclusion": (n_excluded > 0) or (n_flagged > 0),
                "_rank": (
                    target_rank.get(str(target), len(targets)),
                    group_rank.get(str(group), len(group_rank)),
                    sample_rank.get(str(sample), len(sample_rank)),
                ),
            }
        )

    blocks.sort(key=lambda b: (not b["has_exclusion"], b["_rank"]))
    for b in blocks:
        b.pop("_rank", None)
    return blocks


def summarize_dataset(
    df: pd.DataFrame,
    filename: str | None = None,
    analysis_time: datetime | None = None,
) -> dict:
    """Build the upload-time dataset summary the right pane displays.

    Args:
        df: Standardised DataFrame (post :func:`apply_mapping`).
        filename: Original upload filename, for display.
        analysis_time: Override clock (used in tests). Defaults to ``now``.

    Returns:
        dict with the fields the UI renders. ``samples``, ``targets``,
        ``groups``, ``batches`` carry both the count and the list. NaN /
        excluded wells are counted from the standardised data.
    """
    when = (analysis_time or datetime.now()).strftime("%Y-%m-%d %H:%M")
    n_wells = int(len(df))
    samples = list(pd.unique(df["Sample"].astype(str)))
    targets = target_order(df)
    groups = list(pd.unique(df["Group"].astype(str))) if "Group" in df.columns else []

    has_batch = "Batch" in df.columns
    batches = list(pd.unique(df["Batch"].astype(str))) if has_batch else []

    n_na = int(df["Cq"].isna().sum()) if "Cq" in df.columns else 0
    n_excluded = (
        int(df["Excluded"].sum()) if "Excluded" in df.columns else 0
    )

    rep_counts = (
        df.groupby(["Sample", "Target"], sort=False).size().tolist()
        if {"Sample", "Target"}.issubset(df.columns)
        else []
    )
    if rep_counts:
        rep_min, rep_max = min(rep_counts), max(rep_counts)
    else:
        rep_min = rep_max = 0

    return {
        "analysis_time": when,
        "filename": filename,
        "n_wells": n_wells,
        "n_samples": len(samples),
        "samples": samples,
        "n_targets": len(targets),
        "targets": targets,
        "n_groups": len(groups),
        "groups": groups,
        "has_batch_column": has_batch,
        "n_batches": len(batches),
        "batches": batches,
        "n_na_wells": n_na,
        "n_excluded_wells": n_excluded,
        "n_replicate_blocks": len(rep_counts),
        "replicates_min": rep_min,
        "replicates_max": rep_max,
    }
