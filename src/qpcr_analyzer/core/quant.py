from __future__ import annotations

import numpy as np
import pandas as pd


def compute_mean_cq(df: pd.DataFrame) -> pd.DataFrame:
    """Mean Cq per (Target, Group, Sample), dropping wells with Excluded=True or NaN Cq."""
    data = df
    if "Excluded" in data.columns:
        data = data.loc[~data["Excluded"]]
    data = data.dropna(subset=["Cq"])
    grouped = (
        data.groupby(["Target", "Group", "Sample"], sort=False, as_index=False)["Cq"]
        .mean()
        .rename(columns={"Cq": "Mean_Cq"})
    )
    return grouped


def compute_relative_expression(
    df: pd.DataFrame,
    ref_genes: list[str],
    reference_samples: dict[str, list[str]],
) -> dict[str, pd.DataFrame]:
    """Group-aware ΔΔCq relative quantification.

    For each housekeeping gene ``ref`` and each group, relative expression is
    normalized against the mean ΔCq of that group's reference samples, so the
    reference samples in every group anchor at 1.

    Args:
      df: data with columns Well, Target, Sample, Group, Cq and optional Excluded.
      ref_genes: one or more housekeeping genes (names found in ``Target``).
      reference_samples: mapping {group_name: [sample_name, ...]} of the anchor
        samples for each group.

    Returns:
      dict keyed by housekeeping gene → per-target result DataFrame.
    """
    if not ref_genes:
        raise ValueError("At least one housekeeping gene is required.")

    mean_cq = compute_mean_cq(df)
    if mean_cq.empty:
        raise ValueError("No data remains after excluding invalid wells.")
    if mean_cq["Mean_Cq"].isna().any():
        raise ValueError("NaN values present in mean Cq data.")

    available_targets = set(mean_cq["Target"].unique())
    missing = [r for r in ref_genes if r not in available_targets]
    if missing:
        raise ValueError(f"Housekeeping gene(s) not in data: {missing}")

    results: dict[str, pd.DataFrame] = {}
    for ref in ref_genes:
        ref_df = (
            mean_cq.loc[mean_cq["Target"] == ref, ["Sample", "Group", "Mean_Cq"]]
            .rename(columns={"Mean_Cq": f"Mean_Cq_{ref}"})
        )
        target_df = mean_cq.loc[mean_cq["Target"] != ref].copy()
        merged = target_df.merge(ref_df, on=["Sample", "Group"], how="left")
        if merged[f"Mean_Cq_{ref}"].isna().any():
            bad = merged.loc[merged[f"Mean_Cq_{ref}"].isna(), ["Sample", "Group"]]
            raise ValueError(
                f"Housekeeping gene '{ref}' is missing for: "
                f"{bad.drop_duplicates().to_dict(orient='records')}"
            )
        merged["dCq"] = merged["Mean_Cq"] - merged[f"Mean_Cq_{ref}"]
        merged["Is_Reference"] = [
            s in reference_samples.get(g, []) for s, g in zip(merged["Sample"], merged["Group"])
        ]

        anchor = (
            merged.loc[merged["Is_Reference"]]
            .groupby(["Group", "Target"], sort=False, as_index=False)["dCq"]
            .mean()
            .rename(columns={"dCq": "Ref_dCq"})
        )
        merged = merged.merge(anchor, on=["Group", "Target"], how="left")
        missing_groups = (
            merged.loc[merged["Ref_dCq"].isna(), "Group"].drop_duplicates().tolist()
        )
        if missing_groups:
            raise ValueError(
                f"No reference sample designated (or remaining) for group(s): {missing_groups}"
            )
        merged["ddCq"] = merged["dCq"] - merged["Ref_dCq"]
        merged["Relative_Expr"] = np.power(2.0, -merged["ddCq"])
        merged["Reference_Gene"] = ref

        results[ref] = merged[
            [
                "Target",
                "Group",
                "Sample",
                "Mean_Cq",
                f"Mean_Cq_{ref}",
                "dCq",
                "Ref_dCq",
                "ddCq",
                "Relative_Expr",
                "Is_Reference",
                "Reference_Gene",
            ]
        ].reset_index(drop=True)
    return results
