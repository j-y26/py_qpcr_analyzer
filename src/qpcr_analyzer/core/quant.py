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

from .summary import group_order, sample_order, target_order


def _sort_results_by_appearance(
    result_df: pd.DataFrame, source_df: pd.DataFrame
) -> pd.DataFrame:
    """Order result rows by Target → Group → Sample, file-appearance ranks.

    Ranks are derived from ``source_df`` (the standardised input) so the
    ordering matches the order labels first appear in the upload. Rows whose
    Target / Group / Sample is missing from the source — should not happen
    in practice, but defensively handled — are pushed to the end.
    """
    target_rank = {t: i for i, t in enumerate(target_order(source_df))}
    group_rank = {g: i for i, g in enumerate(group_order(source_df))}
    sample_rank = {s: i for i, s in enumerate(sample_order(source_df))}
    out = result_df.copy()
    out["__t_rank"] = out["Target"].astype(str).map(target_rank).fillna(len(target_rank))
    out["__g_rank"] = out["Group"].astype(str).map(group_rank).fillna(len(group_rank))
    out["__s_rank"] = out["Sample"].astype(str).map(sample_rank).fillna(len(sample_rank))
    out = out.sort_values(
        by=["__t_rank", "__g_rank", "__s_rank"], kind="mergesort"
    ).drop(columns=["__t_rank", "__g_rank", "__s_rank"])
    return out.reset_index(drop=True)


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


def samples_missing_hk(
    df: pd.DataFrame, ref_genes: list[str]
) -> dict[str, list[str]]:
    """Return samples that have no usable HK Cq for each housekeeping gene.

    "No usable HK Cq" means: every well of that sample × HK pair is either
    NaN or marked Excluded, so the sample's mean HK Cq cannot be computed.
    The UI uses this to surface affected samples to the user for manual
    confirmation before silently dropping them from analysis.

    Args:
        df: DataFrame with columns Sample, Target, Cq and optional Excluded.
        ref_genes: housekeeping gene names to check.

    Returns:
        dict keyed by HK gene → list of sample names with no valid HK Cq.
        Genes with no missing samples are still present with an empty list.
    """
    mean_cq = compute_mean_cq(df)
    # File-appearance order, not alphanumeric: biological sample labels
    # often encode meaning that lexical sort destroys.
    all_samples = list(pd.unique(df["Sample"].astype(str)))
    out: dict[str, list[str]] = {}
    for ref in ref_genes:
        present = set(
            mean_cq.loc[mean_cq["Target"] == ref, "Sample"].astype(str)
        )
        out[ref] = [s for s in all_samples if s not in present]
    return out


def compute_delta_ct(
    df: pd.DataFrame,
    ref_genes: list[str],
    sample_excludes_per_hk: dict[str, set[str] | list[str]] | None = None,
) -> dict[str, pd.DataFrame]:
    """ΔCt = mean_Cq(Target) − mean_Cq(HK) for each Sample × Target.

    This method normalises expression to the housekeeping gene only; no
    reference group is involved.  The resulting ``dCt`` values can be
    interpreted directly or converted to ``2^(−dCt)`` (``Expr_vs_HK``).

    Args:
        df: DataFrame with columns Well, Target, Sample, Group, Cq and
            optional Excluded.
        ref_genes: housekeeping gene names present in the Target column.
        sample_excludes_per_hk: ``{hk_gene: {sample, ...}}`` — samples to
            drop from each HK's ΔCt sheet (e.g. when a sample lacks a valid
            HK Cq, or the user wants to exclude it for that HK only). The
            sample still appears for other HKs unless listed there too.
            Use the special key ``"*"`` to exclude a sample from every HK.

    Returns:
        dict keyed by housekeeping gene name → DataFrame with columns
        Target, Group, Sample, Mean_Cq, Mean_Cq_{ref}, dCt, Expr_vs_HK,
        Reference_Gene.
    """
    if not ref_genes:
        raise ValueError("At least one housekeeping gene is required.")

    excl_per_hk = {k: set(v) for k, v in (sample_excludes_per_hk or {}).items()}
    global_excl = excl_per_hk.pop("*", set())

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
        drop_samples = global_excl | excl_per_hk.get(ref, set())

        ref_cq = (
            mean_cq.loc[mean_cq["Target"] == ref, ["Sample", "Group", "Mean_Cq"]]
            .rename(columns={"Mean_Cq": ref_col})
        )
        if drop_samples:
            ref_cq = ref_cq[~ref_cq["Sample"].astype(str).isin(drop_samples)]

        target_df = mean_cq.loc[mean_cq["Target"] != ref].copy()
        if drop_samples:
            target_df = target_df[~target_df["Sample"].astype(str).isin(drop_samples)]

        merged = target_df.merge(ref_cq, on=["Sample", "Group"], how="left")
        absent = merged[ref_col].isna()
        if absent.any():
            bad = merged.loc[absent, ["Sample", "Group"]].drop_duplicates()
            raise ValueError(
                f"Housekeeping gene '{ref}' missing for samples: "
                f"{bad.to_dict(orient='records')}. "
                "Add these samples to sample_excludes_per_hk to skip them."
            )
        merged["dCt"] = merged["Mean_Cq"] - merged[ref_col]
        merged["Expr_vs_HK"] = np.power(2.0, -merged["dCt"])
        merged["Reference_Gene"] = ref
        ordered = merged[
            ["Target", "Group", "Sample", "Mean_Cq", ref_col,
             "dCt", "Expr_vs_HK", "Reference_Gene"]
        ]
        results[ref] = _sort_results_by_appearance(ordered, df)
    return results


def compute_delta_delta_ct(
    df: pd.DataFrame,
    ref_genes: list[str],
    reference_group: str,
    sample_batches: dict[str, str] | None = None,
    sample_excludes_per_hk: dict[str, set[str] | list[str]] | None = None,
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
        sample_excludes_per_hk: forwarded to :func:`compute_delta_ct` — see
            its docstring for semantics.

    Returns:
        dict keyed by housekeeping gene name → DataFrame with columns
        Target, Group, Sample, Batch, Mean_Cq, Mean_Cq_{ref},
        dCt, Ref_dCt, ddCt, Relative_Expr,
        Is_Reference_Group, Reference_Gene, Reference_Group.
    """
    if not ref_genes:
        raise ValueError("At least one housekeeping gene is required.")

    dct_results = compute_delta_ct(df, ref_genes, sample_excludes_per_hk)

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
        ordered = merged[
            [
                "Target", "Group", "Sample", "Batch",
                "Mean_Cq", ref_col,
                "dCt", "Ref_dCt", "ddCt", "Relative_Expr",
                "Is_Reference_Group", "Reference_Gene", "Reference_Group",
            ]
        ]
        results[ref] = _sort_results_by_appearance(ordered, df)
    return results
