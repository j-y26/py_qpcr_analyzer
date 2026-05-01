"""NiceGUI front-end for qPCR Analyzer.

This module is a thin presentation layer over :mod:`qpcr_analyzer.core`.
All scientific logic lives in ``core/``; this file only:

1. Renders a five-step :class:`nicegui.ui.stepper`:
   upload → column mapping → groups & batches → reference + outliers → results.
2. Carries per-session state in a plain ``dict`` returned by
   :func:`_new_state` (NiceGUI binds one ``index()`` page per browser
   connection, so this dict is effectively per-session).
3. Calls into ``core`` to do the heavy lifting and renders the returned
   DataFrames as Plotly bar charts plus an Excel download button.

Run with :func:`start` or via the ``qpcr-analyzer`` console script (see
:mod:`qpcr_analyzer.__main__`).
"""

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
    compute_delta_ct,
    compute_delta_delta_ct,
    detect_columns,
    mark_outliers,
    read_table,
    results_to_xlsx_bytes,
    validate_sample_groups,
)
from qpcr_analyzer.core.columns import REQUIRED

PRIMARY = "#2563eb"
PALETTE = ["#2563eb", "#0ea5e9", "#14b8a6", "#f59e0b", "#ef4444", "#a855f7"]


def _new_state() -> dict:
    """Build a fresh per-session state dict.

    Each browser connection gets its own ``state`` instance; nothing is
    shared between users. Keys mirror the five UI steps so each step can
    read what previous steps produced.
    """
    return {
        "filename": None,
        "raw_df": None,
        "mapping": None,
        "standardized": None,
        "sample_batches": {},       # {sample: batch_label}
        "reference_group": None,    # group used as ΔΔCt anchor
        "tolerance": 1.0,
        "flagged": None,
        "excluded_wells": set(),
        "_outliers_initialized": False,
        "ref_genes": [],
        "dct_results": None,
        "ddct_results": None,
    }


@ui.page("/")
def index() -> None:
    """Render the single-page stepper UI.

    NiceGUI calls this function once per browser connection. ``state``
    holds analysis data and ``refs`` holds handles to UI widgets that
    later render-helpers need to mutate (e.g. the outlier table after the
    user changes the tolerance).
    """
    state = _new_state()
    refs: dict[str, Any] = {}

    ui.query("body").classes("bg-slate-50")

    with ui.header(elevated=False).classes(
        "bg-white text-slate-900 items-center border-b border-slate-200 px-6"
    ):
        with ui.row().classes("items-center gap-3"):
            ui.icon("biotech", size="28px").classes("text-blue-600")
            ui.label("qPCR Analyzer").classes("text-xl font-semibold")
            ui.label("· ΔCt & batch-aware ΔΔCt quantification").classes(
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

            # ── Step 1: Upload ────────────────────────────────────────────────
            with ui.step("upload", title="Upload data", icon="upload_file"):
                ui.markdown(
                    "Upload a qPCR results file.  "
                    "Supported formats: **xlsx, xls, csv, tsv, txt**.  "
                    "One row per well."
                ).classes("text-slate-700")

                async def on_upload(e: events.UploadEventArguments) -> None:
                    # NiceGUI 3.x: the event carries a single ``file`` attribute
                    # (FileUpload) with async ``read()`` and a ``name`` field.
                    # Older 2.x exposed ``e.content`` / ``e.name`` directly.
                    try:
                        data = await e.file.read()
                        name = e.file.name
                        df = read_table(io.BytesIO(data), name)
                    except Exception as ex:  # noqa: BLE001
                        ui.notify(f"Failed to read file: {ex}", type="negative")
                        return
                    state["filename"] = name
                    state["raw_df"] = df
                    state["mapping"] = detect_columns(df)
                    refs["file_info"].set_text(
                        f"Loaded {name} — {len(df)} rows × {len(df.columns)} columns"
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

            # ── Step 2: Column mapping ────────────────────────────────────────
            with ui.step("mapping", title="Column mapping", icon="view_column"):
                ui.markdown(
                    "Confirm or adjust which column plays each role.  "
                    "Auto-detected columns show a green badge; required roles are marked."
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
                    # Validate sample↔group consistency
                    sg_errs = validate_sample_groups(std)
                    if sg_errs:
                        for msg in sg_errs:
                            ui.notify(msg, type="negative")
                        return
                    state["standardized"] = std
                    state["_outliers_initialized"] = False
                    state["excluded_wells"] = set()
                    state["sample_batches"] = {}
                    state["reference_group"] = None
                    _render_groups(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Continue", on_click=go_to_groups).props(
                        "color=primary unelevated"
                    )

            # ── Step 3: Groups & batches ──────────────────────────────────────
            with ui.step("groups", title="Groups & batches", icon="groups"):
                ui.markdown(
                    "Verify group assignments and optionally assign samples to "
                    "**batches** (experimental runs).  "
                    "The ΔΔCt method normalises within each batch before merging.  "
                    "Leave all samples in *batch_1* if your data comes from a single run."
                ).classes("text-slate-700")
                refs["groups_container"] = ui.column().classes("w-full gap-3")

                def go_to_outliers() -> None:
                    # Collect batch assignments from the UI widgets
                    batch_map: dict[str, str] = {}
                    for sample, widget in refs.get("batch_inputs", {}).items():
                        val = (widget.value or "").strip() or "batch_1"
                        batch_map[sample] = val
                    state["sample_batches"] = batch_map
                    _render_outliers(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Continue", on_click=go_to_outliers).props(
                        "color=primary unelevated"
                    )

            # ── Step 4: Reference group, HK genes & outliers ─────────────────
            with ui.step("outliers", title="Reference & outlier review", icon="rule"):
                ui.markdown(
                    "Choose the **reference group** (ΔΔCt anchor, expression = 1) "
                    "and **housekeeping gene(s)**, set the replicate tolerance, "
                    "then review flagged wells before running the analysis."
                ).classes("text-slate-700")
                refs["outlier_container"] = ui.column().classes("w-full gap-3")

                def go_to_results() -> None:
                    if state["standardized"] is None:
                        return
                    if not state["ref_genes"]:
                        ui.notify(
                            "Select at least one housekeeping gene.", type="negative"
                        )
                        return
                    ref_group = refs.get("ref_group_sel") and refs["ref_group_sel"].value
                    if not ref_group:
                        ui.notify(
                            "Select a reference group for ΔΔCt.", type="negative"
                        )
                        return
                    state["reference_group"] = ref_group
                    std = state["standardized"].copy()
                    std["Excluded"] = std["Well"].isin(state["excluded_wells"])
                    try:
                        dct = compute_delta_ct(std, list(state["ref_genes"]))
                        ddct = compute_delta_delta_ct(
                            std,
                            list(state["ref_genes"]),
                            reference_group=ref_group,
                            sample_batches=state["sample_batches"] or None,
                        )
                    except Exception as ex:  # noqa: BLE001
                        ui.notify(f"Analysis failed: {ex}", type="negative")
                        return
                    state["dct_results"] = dct
                    state["ddct_results"] = ddct
                    _render_results(state, refs)
                    stepper.next()

                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button("Run analysis", on_click=go_to_results).props(
                        "color=primary unelevated icon=play_arrow"
                    )

            # ── Step 5: Results ───────────────────────────────────────────────
            with ui.step("results", title="Results", icon="bar_chart"):
                refs["results_container"] = ui.column().classes("w-full gap-4")
                with ui.stepper_navigation():
                    ui.button("Back", on_click=stepper.previous).props("flat")
                    ui.button(
                        "Start over",
                        on_click=lambda: ui.navigate.reload(),
                        icon="refresh",
                    ).props("flat")


# ── Render helpers ────────────────────────────────────────────────────────────

def _render_mapping(state: dict, refs: dict) -> None:
    container = refs["mapping_container"]
    container.clear()
    cols = list(state["raw_df"].columns)
    options = ["(none)"] + cols
    mapping = state["mapping"]

    with container:
        for role in ROLES:
            current = mapping.assignments.get(role)
            conf = mapping.confidence.get(role, 0.0)
            required = role in REQUIRED
            with ui.row().classes("w-full items-center gap-3"):
                ui.label(ROLE_LABELS[role]).classes("w-28 font-medium")
                sel = ui.select(
                    options=options,
                    value=current if current else "(none)",
                ).classes("flex-grow").props("outlined dense")

                def _on_change(e, role=role, sel=sel) -> None:  # noqa: ARG001
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
                    ui.badge(f"auto · {conf:.0%}", color="green").classes("text-white")
                elif required:
                    ui.badge("unmatched", color="red").classes("text-white")
                else:
                    ui.badge("n/a", color="grey").classes("text-white")


def _render_groups(state: dict, refs: dict) -> None:
    container = refs["groups_container"]
    container.clear()
    df = state["standardized"]

    samples_meta = (
        df[["Sample", "Group"]]
        .drop_duplicates()
        .sort_values(["Group", "Sample"])
        .reset_index(drop=True)
    )
    groups = samples_meta["Group"].unique().tolist()

    refs["batch_inputs"] = {}

    with container:
        ui.label(
            f"{len(samples_meta)} unique sample(s) across {len(groups)} group(s)"
        ).classes("text-sm text-slate-600")

        with ui.card().classes("w-full border border-slate-200 shadow-none"):
            columns = [
                {"name": "Sample", "label": "Sample", "field": "Sample",
                 "align": "left", "sortable": True},
                {"name": "Group", "label": "Group", "field": "Group",
                 "align": "left", "sortable": True},
                {"name": "Batch", "label": "Batch", "field": "Batch",
                 "align": "left"},
            ]
            rows = [
                {"Sample": r.Sample, "Group": r.Group, "Batch": "batch_1"}
                for r in samples_meta.itertuples()
            ]
            table = ui.table(
                columns=columns, rows=rows, row_key="Sample"
            ).classes("w-full")
            table.add_slot(
                "body",
                r"""
                <q-tr :props="props">
                  <q-td key="Sample" :props="props">{{ props.row.Sample }}</q-td>
                  <q-td key="Group"  :props="props">{{ props.row.Group }}</q-td>
                  <q-td key="Batch"  :props="props">
                    <q-input
                      v-model="props.row.Batch"
                      dense outlined
                      style="min-width:120px"
                      @update:model-value="() => $emit('batch-change', {sample: props.row.Sample, batch: props.row.Batch})"
                    />
                  </q-td>
                </q-tr>
                """,
            )

            # Capture batch edits emitted from the Quasar template
            def _on_batch_change(e) -> None:
                sample = e.args.get("sample")
                batch = (e.args.get("batch") or "").strip() or "batch_1"
                if sample:
                    state["sample_batches"][str(sample)] = batch

            table.on("batch-change", _on_batch_change)

        ui.label(
            "Leave all batch fields as 'batch_1' if your samples come "
            "from a single experimental run."
        ).classes("text-xs text-slate-400 mt-1")


def _render_outliers(state: dict, refs: dict) -> None:
    container = refs["outlier_container"]
    container.clear()
    df = state["standardized"]
    groups = df["Group"].unique().tolist()
    targets = df["Target"].unique().tolist()

    with container:
        with ui.row().classes("w-full items-end gap-4 flex-wrap"):
            ref_group_sel = ui.select(
                options=groups,
                label="Reference group (ΔΔCt anchor)",
                value=state.get("reference_group") or groups[0],
            ).classes("flex-grow min-w-48").props("outlined dense")
            refs["ref_group_sel"] = ref_group_sel

            ref_gene_sel = ui.select(
                options=targets,
                label="Housekeeping gene(s)",
                multiple=True,
                value=state["ref_genes"],
            ).classes("flex-grow min-w-64").props("outlined dense use-chips")

            tol_input = ui.number(
                label="Replicate tolerance (cycles)",
                value=state["tolerance"],
                min=0.1,
                step=0.1,
                format="%.2f",
            ).classes("w-48").props("outlined dense")

            def _apply() -> None:
                try:
                    new_tol = float(tol_input.value)
                except (TypeError, ValueError):
                    new_tol = 1.0
                if new_tol != state["tolerance"]:
                    state["tolerance"] = new_tol
                    state["_outliers_initialized"] = False
                state["ref_genes"] = list(ref_gene_sel.value or [])
                _rebuild_outlier_table(state, refs)

            ui.button("Apply", on_click=_apply).props("color=primary outline")

        refs["outlier_summary_slot"] = ui.label("").classes("text-sm text-slate-600")
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
    auto_wells = flagged.loc[flagged["Outlier"] & flagged["Cq"].notna(), "Well"].tolist()
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
            "Flag": (
                "⚠ outlier" if bool(r.Outlier) and not pd.isna(r.Cq)
                else ("— NA —" if pd.isna(r.Cq) else "ok")
            ),
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
            columns=columns, rows=rows, row_key="Well", pagination=15
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
                data = results_to_xlsx_bytes(
                    std, state["dct_results"], state["ddct_results"]
                )
                out_name = f"qpcr_results_{Path(filename).stem}.xlsx"
                ui.download(data, out_name)

            ui.button("Download xlsx", icon="download", on_click=_download).props(
                "color=primary unelevated"
            )

        ref_group = state["reference_group"]
        batches = sorted(set(state["sample_batches"].values())) if state["sample_batches"] else ["batch_1"]
        ui.label(
            f"Reference group: {ref_group}  ·  "
            f"Batch(es): {', '.join(batches)}"
        ).classes("text-sm text-slate-500")

        with ui.tabs().classes("w-full") as tabs:
            tab_dct = ui.tab("ΔCt", icon="show_chart")
            tab_ddct = ui.tab("ΔΔCt (relative expression)", icon="bar_chart")

        with ui.tab_panels(tabs, value=tab_dct).classes("w-full"):
            with ui.tab_panel(tab_dct):
                ui.markdown(
                    "**ΔCt** = mean Cq(target) − mean Cq(housekeeping gene), "
                    "per sample.  No reference group normalisation.  "
                    "Lower ΔCt = higher expression relative to the housekeeping gene."
                ).classes("text-sm text-slate-600 mb-2")
                _render_result_cards(
                    state["dct_results"], value_col="dCt",
                    y_label="ΔCt", invert_y=False
                )

            with ui.tab_panel(tab_ddct):
                ui.markdown(
                    f"**ΔΔCt** relative expression: "
                    f"2^(−ΔΔCt), normalised so that *{ref_group}* = 1 within each batch."
                ).classes("text-sm text-slate-600 mb-2")
                _render_result_cards(
                    state["ddct_results"], value_col="Relative_Expr",
                    y_label="Relative expression (2^−ΔΔCt)", invert_y=False
                )


def _render_result_cards(
    results: dict[str, pd.DataFrame],
    value_col: str,
    y_label: str,
    invert_y: bool,
) -> None:
    for ref, res_df in results.items():
        with ui.card().classes("w-full border border-slate-200 shadow-none"):
            with ui.row().classes("items-center gap-3"):
                ui.label(f"Housekeeping: {ref}").classes("text-base font-semibold")
                ui.label(
                    f"{res_df['Target'].nunique()} target(s) · "
                    f"{res_df['Group'].nunique()} group(s) · "
                    f"{res_df['Sample'].nunique()} sample(s)"
                ).classes("text-xs text-slate-500")

            targets = res_df["Target"].unique().tolist()
            groups = res_df["Group"].unique().tolist()
            use_group = len(groups) > 1
            cols = 1 if len(targets) == 1 else 2
            with ui.grid(columns=cols).classes("w-full gap-3"):
                for target in targets:
                    sub = res_df[res_df["Target"] == target]
                    fig = _plot_figure(
                        sub, target, ref, use_group, value_col, y_label, invert_y
                    )
                    ui.plotly(fig).classes("w-full h-80")


def _plot_figure(
    df: pd.DataFrame,
    target: str,
    ref: str,
    use_group: bool,
    value_col: str,
    y_label: str,
    invert_y: bool,
) -> go.Figure:
    x_col = "Group" if use_group else "Sample"
    summary = (
        df.groupby(x_col, sort=False)[value_col]
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
            y=df[value_col],
            mode="markers",
            marker=dict(color="#0f172a", size=6, opacity=0.75),
            name="samples",
            hovertemplate="%{x}<br>%{customdata}<br>val = %{y:.3f}<extra></extra>",
            customdata=df["Sample"],
        )
    fig.update_layout(
        title=dict(text=f"{target} / {ref}", x=0.5, xanchor="center"),
        yaxis_title=y_label,
        yaxis_autorange="reversed" if invert_y else True,
        xaxis_title="",
        template="plotly_white",
        showlegend=False,
        margin=dict(l=40, r=10, t=50, b=40),
        plot_bgcolor="white",
    )
    return fig


def start(host: str | None = None, port: int | None = None) -> None:
    """Launch the NiceGUI server (blocking).

    Args:
        host: Bind address. Falls back to ``$QPCR_HOST`` and finally
            ``"127.0.0.1"`` (loopback only — change to ``"0.0.0.0"`` to
            expose on the LAN).
        port: TCP port. Falls back to ``$QPCR_PORT`` and finally ``8090``.
            Port 8080 is intentionally avoided because Thermo Fisher's
            *Design & Analysis* software (commonly used to inspect qPCR
            raw data on the same machine) listens there.

    The call blocks until the server is stopped. ``reload=False`` is set
    so the package works correctly when installed (NiceGUI's hot-reload
    requires a writable source tree). ``show=False`` prevents NiceGUI from
    auto-opening a browser tab — the user opens the printed URL manually.
    """
    host = host or os.environ.get("QPCR_HOST", "127.0.0.1")
    port = int(port or os.environ.get("QPCR_PORT", 8090))
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
