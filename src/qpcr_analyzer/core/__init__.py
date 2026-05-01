"""Pure-Python analysis library for qPCR data.

The submodules below have **no UI dependency** and can be used standalone:

- :mod:`.io`        file readers (xlsx / xls / csv / tsv / txt → DataFrame).
- :mod:`.columns`   fuzzy column-role detection, mapping, sample/group validation.
- :mod:`.outliers`  replicate outlier flagging using the tightest-cluster rule.
- :mod:`.quant`     ΔCt and batch-aware ΔΔCt relative quantification.
- :mod:`.export`    Excel-workbook writer (raw + per-HK ΔCt/ΔΔCt + Prism sheets).

A typical pipeline::

    df    = read_table(path, path)
    mp    = detect_columns(df)
    std   = apply_mapping(df, mp)
    flags = mark_outliers(std, tolerance=1.0)
    std["Excluded"] = flags["Outlier"]
    dct   = compute_delta_ct(std, ref_genes=["GAPDH"])
    ddct  = compute_delta_delta_ct(std, ["GAPDH"], reference_group="ctrl")
    xlsx  = results_to_xlsx_bytes(std, dct, ddct)
"""

from .columns import (
    ROLE_LABELS,
    ROLES,
    ColumnMapping,
    apply_mapping,
    detect_columns,
    validate_sample_groups,
)
from .export import results_to_xlsx_bytes
from .io import read_table
from .outliers import mark_outliers
from .quant import compute_delta_ct, compute_delta_delta_ct, compute_mean_cq

__all__ = [
    "ROLES",
    "ROLE_LABELS",
    "ColumnMapping",
    "apply_mapping",
    "compute_delta_ct",
    "compute_delta_delta_ct",
    "compute_mean_cq",
    "detect_columns",
    "mark_outliers",
    "read_table",
    "results_to_xlsx_bytes",
    "validate_sample_groups",
]
