from __future__ import annotations

import io
import os
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go
from nicegui import events, ui

from qpcr_analyzer.core import (
    ROLE_LABELS,
    ROLES,
    apply_mapping,
    compute_relative_expression,
    detect_columns,
    mark_outliers,
    read_table,
    results_to_xlsx_bytes,
)
from qpcr_analyzer.core.columns import REQUIRED

PRIMARY = "#2563eb"
PALETTE = ["#2563eb", "#0ea5e9", "#14b8a6", "#f59e0b", "#ef4444", "#a855f7"]


def _new_state() -> dict:
    return {
        "filename": None,
        "raw_df": None,
        "mapping": None,
        "standardized": None,
        "tolerance": 1.0,
        "flagged": None,
        "excluded_wells": set(),
        "_outliers_initialized": False,
        "reference_samples": {},
        "ref_genes": [],
        "results": None,
    }


@ui.page("/")
def index() -> None:
    state = _new_state()
    refs: dict[str, Any] = {}

    ui.query("body").classes("bg-slate-50")

    with ui.header(elevated=False).classes(
        "bg-white text-slate-900 items-center border-b border-slate-200 px-6"
    ):
        with ui.row().classes("items-center gap-3"):
            ui.icon("biotech", size="28px").classes("text-blue-600")
            ui.label("qPCR Analyzer").classes("text-xl font-semibold")
            ui.label("· lightweight relative quantification").classes(
                "text-sm text-slate-500"
            )
        ui.space()
        dark = ui.dark_mode()
        ui.button(icon="dark_mode", on_click=dark.toggle).props("flat round dense")

    with ui.column().classes("w-full max-w-5xl mx-auto p-6 gap-4"):
        with ui.stepper().props("vertical flat header-nav").classes(
            "w-full bg-white rounded-lg shadow-sm border border-slate-200"
        ) as stepper:
            refs["stepper"] = stepper

            with ui.step("upload", title="Upload data", icon="upload_file"):
                ui.markdown(
                    "Upload a qPCR dataset. Supported formats: **xlsx, xls, csv, "
                    "tsv, txt**. One row per well."
                ).classes("text-slate-700")

                async def on_upload(e: events.UploadEventArguments) -> None:
                    try:
                        raw_bytes = e.content.read()
                        df = read_table(io.BytesIO(raw_bytes), e.name)
                    except Exception as ex:  # noqa: BLE001
                        ui.notify(f"Failed to read file: {ex}", type="negative")
                        return
                    state["filename"] = e.name
                    state["raw_df"] = df
                    state["mapping"] = detect_columns(df)
                    refs["file_info"].set_text(
                        f"Loaded {e.name} — {len(df)} rows × {len(df.columns)} columns"
                    )
                    _render_mapping(state, refs)
                    refs["upload_next"].enable()

                ui.upload(on_upload=on_upload, auto_upload=True, max_files=1).props(
                    "accept=.xlsx,.xls,.csv,.tsv,.txt color=primary flat bordered"
                ).classes("w-full")
                refs["file_info"] = ui.label("").classes("text-sm text-slate-600")

                with ui.stepper_navigation():
                    refs["upload_next"] = (
                        ui.button("Continue", on_click=stepper.next)
                        .props("color=primary unelevated")
                    )
                    refs["upload_next"].disable()

            with ui.step("mapping", title="Column mapping", icon="view_column"):
                ui.markdown(
                    "Confirm or adjust how each column is used. Auto-detected "
                    "columns show a green check; required roles are marked."
                ).classes("text-slate-700")
                refs["mapping_container"] = ui.column().classes("w-full gap-2")

                def go_to_groups() -> None:
                    mapping = state["mapping"]
                    if mapping is None:
                        return
                    errs = mapping.validate()
                    if errs:
                        for msg in errs:
                            ui.notify(msg, type="negative")
                        return
                    std = apply_mapping(state["raw_df"], mapping)
                    if std["Cq"].notna().sum() == 0:
                        ui.notify(
                            "Selected Cq column has no numeric values.",
                            type="negative",
                        )
                        return
                    state["standardized"] = std
                    state["_outliers_initialized"] = False
                    state["excluded_wells"] = set()
                    state["reference_samples"] = {}
                    _render_groups(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Continue", on_click=go_to_groups).props(
                        "color=primary unelevated"
                    )

            with ui.step("groups", title="Groups & reference samples", icon="groups"):
                ui.markdown(
                    "Mark the reference samples in each group. Relative expression "
                    "within a group is anchored so the mean of its reference "
                    "samples equals 1. Single-group datasets default to all "
                    "samples as reference."
                ).classes("text-slate-700")
                refs["groups_container"] = ui.column().classes("w-full gap-3")

                def go_to_outliers() -> None:
                    ref_map: dict[str, list[str]] = {}
                    switches = refs.get("group_switches", {})
                    for group, sw_map in switches.items():
                        chosen = [s for s, sw in sw_map.items() if sw.value]
                        if not chosen:
                            ui.notify(
                                f"Group '{group}' has no reference sample "
                                "selected — mark at least one.",
                                type="negative",
                            )
                            return
                        ref_map[group] = chosen
                    state["reference_samples"] = ref_map
                    _render_outliers(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Continue", on_click=go_to_outliers).props(
                        "color=primary unelevated"
                    )

            with ui.step("outliers", title="Outlier review", icon="rule"):
                ui.markdown(
                    "Wells with NaN Cq or replicate disagreement beyond the "
                    "tolerance are pre-selected for exclusion. Adjust as needed, "
                    "then choose housekeeping gene(s) and run the analysis."
                ).classes("text-slate-700")
                refs["outlier_container"] = ui.column().classes("w-full gap-3")

                def go_to_results() -> None:
                    if state["standardized"] is None:
                        return
                    if not state["ref_genes"]:
                        ui.notify(
                            "Select at least one housekeeping gene.",
                            type="negative",
                        )
                        return
                    std = state["standardized"].copy()
                    std["Excluded"] = std["Well"].isin(state["excluded_wells"])
                    try:
                        results = compute_relative_expression(
                            std,
                            list(state["ref_genes"]),
                            state["reference_samples"],
                        )
                    except Exception as ex:  # noqa: BLE001
                        ui.notify(f"Analysis failed: {ex}", type="negative")
                        return
                    state["results"] = results
                    _render_results(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Run analysis", on_click=go_to_results).props(
                        "color=primary unelevated icon=play_arrow"
                    )

            with ui.step("results", title="Results", icon="bar_chart"):
                refs["results_container"] = ui.column().classes("w-full gap-4")
                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button(
                        "Start over",
                        on_click=lambda: ui.navigate.reload(),
                        icon="refresh",
                    ).props("flat")


def _render_mapping(state: dict, refs: dict) -> None:
    container = refs["mapping_container"]
    container.clear()
    cols = list(state["raw_df"].columns)
    options = ["(none)"] + list(cols)
    mapping = state["mapping"]

    with container:
        for role in ROLES:
            current = mapping.assignments.get(role)
            conf = mapping.confidence.get(role, 0.0)
            required = role in REQUIRED
            with ui.row().classes("w-full items-center gap-3 q-mb-xs"):
                ui.label(ROLE_LABELS[role]).classes("w-28 font-medium")
                sel = ui.select(
                    options=options,
                    value=current if current else "(none)",
                ).classes("flex-grow").props("outlined dense")

                def _on_change(e, role=role, sel=sel) -> None:
                    v = sel.value
                    mapping.assignments[role] = None if v == "(none)" else v

                sel.on_value_change(_on_change)

                if required:
                    ui.label("required").classes(
                        "text-xs text-slate-500 uppercase tracking-wide"
                    )
                else:
                    ui.label("optional").classes(
                        "text-xs text-slate-400 uppercase tracking-wide"
                    )

                if current is not None:
                    badge = ui.badge(f"auto · {conf:.0%}", color="green")
                    badge.classes("text-white")
                elif required:
                    ui.badge("unmatched", color="red").classes("text-white")
                else:
                    ui.badge("n/a", color="grey").classes("text-white")


def _render_groups(state: dict, refs: dict) -> None:
    container = refs["groups_container"]
    container.clear()
    df = state["standardized"]
    samples_by_group = (
        df[["Group", "Sample"]]
        .drop_duplicates()
        .groupby("Group", sort=False)["Sample"]
        .apply(list)
        .to_dict()
    )
    single_group = len(samples_by_group) == 1
    refs["group_switches"] = {}

    with container:
        if single_group:
            ui.label(
                "Single group detected — by default every sample is a reference "
                "(relative expression centers on 1). Toggle below if a subset "
                "should anchor instead."
            ).classes("text-sm text-slate-600")

        for group, samples in samples_by_group.items():
            with ui.card().classes("w-full border border-slate-200 shadow-none"):
                with ui.row().classes("items-center w-full"):
                    ui.label(f"Group: {group}").classes("text-base font-semibold")
                    ui.label(f"{len(samples)} sample(s)").classes(
                        "text-xs text-slate-500"
                    )
                    ui.space()

                    def _mark_all(g=group) -> None:
                        for sw in refs["group_switches"][g].values():
                            sw.value = True

                    def _clear_all(g=group) -> None:
                        for sw in refs["group_switches"][g].values():
                            sw.value = False

                    ui.button("All", on_click=_mark_all).props("flat dense size=sm")
                    ui.button("None", on_click=_clear_all).props(
                        "flat dense size=sm"
                    )

                sw_map: dict[str, Any] = {}
                with ui.row().classes("flex-wrap gap-x-4 gap-y-1"):
                    for s in samples:
                        sw = ui.switch(s, value=single_group).props("dense")
                        sw_map[s] = sw
                refs["group_switches"][group] = sw_map


def _render_outliers(state: dict, refs: dict) -> None:
    container = refs["outlier_container"]
    container.clear()
    df = state["standardized"]
    targets = df["Target"].unique().tolist()

    with container:
        with ui.row().classes("w-full items-end gap-4 flex-wrap"):
            tol_input = ui.number(
                label="Replicate tolerance (cycles)",
                value=state["tolerance"],
                min=0.1,
                step=0.1,
                format="%.2f",
            ).classes("w-48").props("outlined dense")
            ref_sel = ui.select(
                options=targets,
                label="Housekeeping gene(s)",
                multiple=True,
                value=state["ref_genes"],
            ).classes("flex-grow min-w-64").props("outlined dense use-chips")

            def _apply() -> None:
                try:
                    new_tol = float(tol_input.value)
                except (TypeError, ValueError):
                    new_tol = 1.0
                if new_tol != state["tolerance"]:
                    state["tolerance"] = new_tol
                    state["_outliers_initialized"] = False
                state["ref_genes"] = list(ref_sel.value or [])
                _rebuild_outlier_table(state, refs)

            ui.button("Apply", on_click=_apply).props("color=primary outline")

        refs["outlier_summary_slot"] = ui.label("").classes(
            "text-sm text-slate-600"
        )
        refs["outlier_table_slot"] = ui.column().classes("w-full")
        refs["excluded_select_slot"] = ui.column().classes("w-full")
        _rebuild_outlier_table(state, refs)


def _rebuild_outlier_table(state: dict, refs: dict) -> None:
    slot = refs["outlier_table_slot"]
    slot.clear()
    sel_slot = refs["excluded_select_slot"]
    sel_slot.clear()

    flagged = mark_outliers(state["standardized"], tolerance=state["tolerance"])
    state["flagged"] = flagged
    auto_wells = flagged.loc[
        flagged["Outlier"] & flagged["Cq"].notna(), "Well"
    ].tolist()
    nan_wells = flagged.loc[flagged["Cq"].isna(), "Well"].tolist()
    if not state.get("_outliers_initialized"):
        state["excluded_wells"] = set(auto_wells) | set(nan_wells)
        state["_outliers_initialized"] = True

    df = flagged.sort_values(
        by=["Outlier", "Target", "Sample", "Well"],
        ascending=[False, True, True, True],
    ).reset_index(drop=True)

    rows = [
        {
            "Well": r.Well,
            "Target": r.Target,
            "Group": r.Group,
            "Sample": r.Sample,
            "Cq": None if pd.isna(r.Cq) else round(float(r.Cq), 3),
            "Replicates": int(r.Replicates),
            "Flag": "⚠ outlier"
            if bool(r.Outlier) and not pd.isna(r.Cq)
            else ("— NA —" if pd.isna(r.Cq) else "ok"),
            "_flagged": bool(r.Outlier),
        }
        for r in df.itertuples()
    ]
    columns = [
        {"name": c, "label": c, "field": c, "align": "left", "sortable": True}
        for c in ["Well", "Target", "Group", "Sample", "Cq", "Replicates", "Flag"]
    ]

    with slot:
        table = ui.table(
            columns=columns,
            rows=rows,
            row_key="Well",
            pagination=15,
        ).classes("w-full")
        table.add_slot(
            "body",
            r"""
            <q-tr :props="props" :class="props.row._flagged ? 'bg-red-1' : ''">
              <q-td v-for="col in props.cols" :key="col.name" :props="props">
                {{ col.value }}
              </q-td>
            </q-tr>
            """,
        )

    all_wells = flagged["Well"].tolist()
    with sel_slot:
        excluded_select = (
            ui.select(
                options=all_wells,
                value=sorted(state["excluded_wells"]),
                multiple=True,
                label="Wells to exclude from analysis",
            )
            .classes("w-full")
            .props("outlined dense use-chips clearable")
        )

        def _on_change(_e) -> None:
            state["excluded_wells"] = set(excluded_select.value or [])
            refs["outlier_summary_slot"].set_text(_excluded_summary(state))

        excluded_select.on_value_change(_on_change)

    refs["outlier_summary_slot"].set_text(_excluded_summary(state))


def _excluded_summary(state: dict) -> str:
    n_total = 0 if state["flagged"] is None else len(state["flagged"])
    n_excl = len(state["excluded_wells"])
    return f"{n_excl} of {n_total} well(s) marked for exclusion."


def _render_results(state: dict, refs: dict) -> None:
    container = refs["results_container"]
    container.clear()
    filename = state["filename"] or "data"

    with container:
        with ui.row().classes("w-full items-center"):
            ui.label("Analysis results").classes("text-lg font-semibold")
            ui.space()

            def _download() -> None:
                std = state["standardized"].copy()
                std["Excluded"] = std["Well"].isin(state["excluded_wells"])
                data = results_to_xlsx_bytes(std, state["results"])
                out_name = f"qpcr_results_{Path(filename).stem}.xlsx"
                ui.download(data, out_name)

            ui.button(
                "Download xlsx", icon="download", on_click=_download
            ).props("color=primary unelevated")

        for ref, res_df in state["results"].items():
            with ui.card().classes("w-full border border-slate-200 shadow-none"):
                with ui.row().classes("items-center gap-3"):
                    ui.label(f"Housekeeping: {ref}").classes(
                        "text-base font-semibold"
                    )
                    ui.label(
                        f"{res_df['Target'].nunique()} target(s) · "
                        f"{res_df['Group'].nunique()} group(s) · "
                        f"{res_df['Sample'].nunique()} sample(s)"
                    ).classes("text-xs text-slate-500")

                groups = res_df["Group"].unique().tolist()
                use_group = len(groups) > 1
                targets = res_df["Target"].unique().tolist()
                cols = 1 if len(targets) == 1 else 2
                with ui.grid(columns=cols).classes("w-full gap-3"):
                    for target in targets:
                        sub = res_df[res_df["Target"] == target]
                        fig = _plot_figure(sub, target, ref, use_group)
                        ui.plotly(fig).classes("w-full h-80")


def _plot_figure(
    df: pd.DataFrame,
    target: str,
    ref: str,
    use_group: bool,
) -> go.Figure:
    x_col = "Group" if use_group else "Sample"
    summary = (
        df.groupby(x_col, sort=False)["Relative_Expr"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .fillna(0)
    )
    fig = go.Figure()
    fig.add_bar(
        x=summary[x_col],
        y=summary["mean"],
        error_y=dict(
            type="data", array=summary["std"], visible=True, color="#334155"
        ),
        marker=dict(color=PRIMARY, line=dict(color="#1e3a8a", width=0)),
        name="mean",
        hovertemplate=f"%{{x}}<br>mean = %{{y:.3f}}<extra></extra>",
    )
    if use_group:
        fig.add_scatter(
            x=df[x_col],
            y=df["Relative_Expr"],
            mode="markers",
            marker=dict(color="#0f172a", size=6, opacity=0.75),
            name="samples",
            hovertemplate="%{x}<br>%{customdata}<br>rel = %{y:.3f}<extra></extra>",
            customdata=df["Sample"],
        )
    fig.update_layout(
        title=dict(text=f"{target} / {ref}", x=0.5, xanchor="center"),
        yaxis_title=f"{target} / {ref}",
        xaxis_title="",
        template="plotly_white",
        showlegend=False,
        margin=dict(l=40, r=10, t=50, b=40),
        plot_bgcolor="white",
    )
    return fig


def start(host: str | None = None, port: int | None = None) -> None:
    host = host or os.environ.get("QPCR_HOST", "127.0.0.1")
    port = int(port or os.environ.get("QPCR_PORT", 8080))
    ui.run(
        title="qPCR Analyzer",
        host=host,
        port=port,
        reload=False,
        show=False,
        favicon="🧬",
    )


if __name__ in {"__main__", "__mp_main__"}:
    start()
