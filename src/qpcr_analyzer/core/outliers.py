"""Replicate-outlier flagging.

For every (Sample, Target) group, find the tightest cluster of replicates
whose Cq range is within a user tolerance (default 1 cycle) and flag the
rest. NaN Cq values are always flagged. The cluster rule is a generalisation
of the classic R "closest-pair" triplicate rule that works for any number
of replicates ≥ 1.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _flag_cluster(values: np.ndarray, tol: float) -> np.ndarray:
    """Find the tightest contiguous-in-sorted-order cluster whose range <= tol.

    Generalizes the original R "closest-pair" triplicate rule to any n >= 1:
      - n == 1: flag the well (no replicate to confirm against).
      - n == 2: keep both if |a - b| <= tol, else flag both.
      - n >= 3: keep the longest sorted run whose range <= tol; flag the rest.
        Ties on length are broken by tightest range. A true tie (same length
        and same range) is treated as ambiguous and flags all wells.
      - If the longest valid run has fewer than 2 members, flag all.
    """
    n = len(values)
    flags = np.ones(n, dtype=bool)
    if n <= 1:
        return flags

    order = np.argsort(values, kind="mergesort")
    sorted_vals = values[order]

    best_len = 1
    best_starts: list[int] = []
    for i in range(n):
        j = i
        while j + 1 < n and sorted_vals[j + 1] - sorted_vals[i] <= tol:
            j += 1
        run_len = j - i + 1
        if run_len > best_len:
            best_len = run_len
            best_starts = [i]
        elif run_len == best_len:
            best_starts.append(i)

    if best_len < 2:
        return flags

    ranges = [
        (sorted_vals[s + best_len - 1] - sorted_vals[s], s) for s in best_starts
    ]
    ranges.sort()
    if len(ranges) > 1 and ranges[0][0] == ranges[1][0]:
        return flags

    start = ranges[0][1]
    clean_idx = order[start : start + best_len]
    flags[clean_idx] = False
    return flags


def mark_outliers(df: pd.DataFrame, tolerance: float = 1.0) -> pd.DataFrame:
    """Add Outlier (bool) and Replicates (int) columns.

    Group rows by (Sample, Target). NaN Cq is always flagged. Remaining Cq
    values are evaluated with the cluster rule in ``_flag_cluster``.
    """
    out = df.copy()
    out["Outlier"] = False
    out["Replicates"] = 0

    groups = out.groupby(["Sample", "Target"], sort=False).indices
    for _, idx in groups.items():
        cq = out["Cq"].to_numpy(dtype=float)[idx]
        nan_mask = np.isnan(cq)
        flags = np.ones(len(cq), dtype=bool)
        if (~nan_mask).any():
            valid_flags = _flag_cluster(cq[~nan_mask], tolerance)
            flags[~nan_mask] = valid_flags
        out.iloc[idx, out.columns.get_loc("Outlier")] = flags
        out.iloc[idx, out.columns.get_loc("Replicates")] = int((~nan_mask).sum())

    return out
