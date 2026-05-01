"""ΔCt and batch-aware ΔΔCt relative quantification.

The two top-level functions both consume a *standardised* DataFrame produced
by :func:`qpcr_analyzer.core.columns.apply_mapping` (columns ``Well``,
``Target``, ``Sample``, ``Group``, ``Cq``, optional ``Excluded``) and return
a ``dict`` keyed by housekeeping-gene name → result DataFrame, so callers
can iterate through one sheet per HK gene without re-running the pipeline.

Why two functions?
    * :func:`compute_delta_ct` normalises against the housekeeping gene
      only; no biological reference group is required. Useful for datasets
      without a clear control or for sample-by-sample HK-normalised values.
    * :func:`compute_delta_delta_ct` additionally normalises against a
      *reference biological group*, anchoring its mean ΔCt at relative
      expression = 1 within each batch. The exact invariant is
      ``mean(ΔΔCt_ref) = 0`` per batch — ``mean(2^(−ΔΔCt)_ref) = 1`` only
      when the reference samples have identical ΔCt, because exponentials
      do not commute with averaging.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_mean_cq(df: pd.DataFrame) -> pd.DataFrame:
    """Mean Cq per (Target, Group, Sample), dropping excluded/NaN wells."""
    data = df
    if "Excluded" in data.columns:
        data = data.loc[~data["Excluded"]]
    data = data.dropna(subset=["Cq"])
    return (
        data.groupby(["Target", "Group", "Sample"], sort=False, as_index=False)["Cq"]
        .mean()
        .rename(columns={"Cq": "Mean_Cq"})
    )


def compute_delta_ct(
    df: pd.DataFrame,
    ref_genes: list[str],
) -> dict[str, pd.DataFrame]:
    """ΔCt = mean_Cq(Target) − mean_Cq(HK) for each Sample × Target.

    This method normalises expression to the housekeeping gene only; no
    reference group is involved.  The resulting ``dCt`` values can be
    interpreted directly or converted to ``2^(−dCt)`` (``Expr_vs_HK``).

    Args:
        df: DataFrame with columns Well, Target, Sample, Group, Cq and
            optional Excluded.
        ref_genes: housekeeping gene names present in the Target column.

    Returns:
        dict keyed by housekeeping gene name → DataFrame with columns
        Target, Group, Sample, Mean_Cq, Mean_Cq_{ref}, dCt, Expr_vs_HK,
        Reference_Gene.
    """
    if not ref_genes:
        raise ValueError("At least one housekeeping gene is required.")

    mean_cq = compute_mean_cq(df)
    if mean_cq.empty:
        raise ValueError("No data remains after excluding invalid wells.")

    available = set(mean_cq["Target"].unique())
    missing = [r for r in ref_genes if r not in available]
    if missing:
        raise ValueError(f"Housekeeping gene(s) not found in data: {missing}")

    results: dict[str, pd.DataFrame] = {}
    for ref in ref_genes:
        ref_col = f"Mean_Cq_{ref}"
        ref_cq = (
            mean_cq.loc[mean_cq["Target"] == ref, ["Sample", "Group", "Mean_Cq"]]
            .rename(columns={"Mean_Cq": ref_col})
        )
        target_df = mean_cq.loc[mean_cq["Target"] != ref].copy()
        merged = target_df.merge(ref_cq, on=["Sample", "Group"], how="left")
        absent = merged[ref_col].isna()
        if absent.any():
            bad = merged.loc[absent, ["Sample", "Group"]].drop_duplicates()
            raise ValueError(
                f"Housekeeping gene '{ref}' missing for samples: "
                f"{bad.to_dict(orient='records')}"
            )
        merged["dCt"] = merged["Mean_Cq"] - merged[ref_col]
        merged["Expr_vs_HK"] = np.power(2.0, -merged["dCt"])
        merged["Reference_Gene"] = ref
        results[ref] = merged[
            ["Target", "Group", "Sample", "Mean_Cq", ref_col,
             "dCt", "Expr_vs_HK", "Reference_Gene"]
        ].reset_index(drop=True)
    return results


def compute_delta_delta_ct(
    df: pd.DataFrame,
    ref_genes: list[str],
    reference_group: str,
    sample_batches: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Batch-aware ΔΔCt relative quantification.

    Within each batch the mean ΔCt of *reference_group* samples anchors at
    relative expression = 1.  Results from all batches are merged so that
    cross-batch comparisons can be visualised together.

    Args:
        df: DataFrame with columns Well, Target, Sample, Group, Cq and
            optional Excluded.
        ref_genes: housekeeping gene names present in the Target column.
        reference_group: the Group value whose samples serve as the ΔΔCt
            denominator (their batch-mean relative expression = 1).
        sample_batches: mapping {sample_name: batch_label}.  When None or a
            sample is absent from the dict, that sample is assigned to
            "batch_1" (i.e. all samples belong to a single batch by default).

    Returns:
        dict keyed by housekeeping gene name → DataFrame with columns
        Target, Group, Sample, Batch, Mean_Cq, Mean_Cq_{ref},
        dCt, Ref_dCt, ddCt, Relative_Expr,
        Is_Reference_Group, Reference_Gene, Reference_Group.
    """
    if not ref_genes:
        raise ValueError("At least one housekeeping gene is required.")

    dct_results = compute_delta_ct(df, ref_genes)

    all_samples = df["Sample"].unique()
    if sample_batches is None:
        batch_map: dict[str, str] = {str(s): "batch_1" for s in all_samples}
    else:
        batch_map = {str(s): sample_batches.get(str(s), "batch_1") for s in all_samples}

    results: dict[str, pd.DataFrame] = {}
    for ref, dct_df in dct_results.items():
        dct_df = dct_df.copy()
        dct_df["Batch"] = dct_df["Sample"].map(batch_map).fillna("batch_1")
        dct_df["Is_Reference_Group"] = dct_df["Group"] == reference_group

        ref_rows = dct_df[dct_df["Is_Reference_Group"]]
        if ref_rows.empty:
            available_groups = dct_df["Group"].unique().tolist()
            raise ValueError(
                f"Reference group '{reference_group}' not found in data. "
                f"Available groups: {available_groups}"
            )

        # Per (Batch × Target): mean ΔCt of the reference group
        anchor = (
            ref_rows.groupby(["Batch", "Target"], sort=False, as_index=False)["dCt"]
            .mean()
            .rename(columns={"dCt": "Ref_dCt"})
        )

        # Verify every (Batch × Target) combination has a reference anchor
        all_combos = dct_df[["Batch", "Target"]].drop_duplicates()
        check = all_combos.merge(anchor, on=["Batch", "Target"], how="left")
        missing = check[check["Ref_dCt"].isna()]
        if not missing.empty:
            pairs = missing[["Batch", "Target"]].to_dict(orient="records")
            raise ValueError(
                f"Reference group '{reference_group}' has no samples in "
                f"the following (batch, target) combinations: {pairs}. "
                "Ensure every batch contains at least one reference-group sample "
                "for every measured target."
            )

        merged = dct_df.merge(anchor, on=["Batch", "Target"], how="left")
        merged["ddCt"] = merged["dCt"] - merged["Ref_dCt"]
        merged["Relative_Expr"] = np.power(2.0, -merged["ddCt"])
        merged["Reference_Gene"] = ref
        merged["Reference_Group"] = reference_group

        ref_col = f"Mean_Cq_{ref}"
        results[ref] = merged[
            [
                "Target", "Group", "Sample", "Batch",
                "Mean_Cq", ref_col,
                "dCt", "Ref_dCt", "ddCt", "Relative_Expr",
                "Is_Reference_Group", "Reference_Gene", "Reference_Group",
            ]
        ].reset_index(drop=True)
    return results
