"""NiceGUI front-end for qPCR Analyzer.

The page is a two-pane layout:

* **Left pane** — a vertical stepper that drives the workflow:
  upload → column mapping → groups & batches → outlier review →
  housekeeping gene selection → run ΔCt → reference group & run ΔΔCt.
* **Right pane** — persistent output tabs that activate as data flows in:
  *Summary*, *Data preview*, *Excluded blocks*, *ΔCt results*, *ΔΔCt
  results*, *Downloads*.

All scientific logic lives in :mod:`qpcr_analyzer.core`; this module only
collects user input, calls the core, and renders Plotly figures + tables.
"""

from __future__ import annotations

import asyncio
import io
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go
from nicegui import events, ui
from plotly.subplots import make_subplots

from qpcr_analyzer import __version__
from qpcr_analyzer.core import (
    ROLE_LABELS,
    ROLES,
    apply_mapping,
    build_blocks,
    compute_delta_ct,
    compute_delta_delta_ct,
    detect_columns,
    group_order,
    mark_outliers,
    read_table,
    results_to_csv_zip_bytes,
    results_to_xlsx_bytes,
    sample_order,
    samples_missing_hk,
    sort_wells,
    summarize_dataset,
    target_order,
    validate_sample_batches,
    validate_sample_groups,
)
from qpcr_analyzer.core.columns import REQUIRED

PRIMARY = "#2563eb"
# Distinct, colour-blind-friendly palette used for per-bar colouring.
PALETTE = [
    "#2563eb", "#ef4444", "#10b981", "#f59e0b", "#a855f7", "#0ea5e9",
    "#ec4899", "#14b8a6", "#84cc16", "#f97316", "#6366f1", "#dc2626",
    "#06b6d4", "#eab308", "#8b5cf6", "#22c55e",
]


def _build_color_map(items: list[str]) -> dict[str, str]:
    """Stable colour-by-label mapping that wraps the palette if needed."""
    return {str(it): PALETTE[i % len(PALETTE)] for i, it in enumerate(items)}


def _filename_stem(state: dict) -> str:
    """`qpcr_analysis_<YYYYMMDD>_<original-stem>` — used for downloads."""
    stem = Path(state.get("filename") or "data").stem
    today = datetime.now().strftime("%Y%m%d")
    return f"qpcr_analysis_{today}_{stem}"


def _new_state() -> dict:
    """Build a fresh per-session state dict."""
    return {
        "filename": None,
        "raw_df": None,
        "mapping": None,
        "standardized": None,
        "summary": None,
        "sample_batches": {},        # {sample: batch_label}
        "has_batches": False,        # toggled by user via the step-3 checkbox
        "reference_group": None,
        "tolerance": 1.0,
        "flagged": None,
        "excluded_wells": set(),
        "_outliers_initialized": False,
        "ref_genes": [],
        # {hk_gene: set(sample_name, ...)}, "*" means all housekeeping genes
        "sample_excludes_per_hk": {},
        "hk_applied": False,
        "dct_results": None,
        "ddct_results": None,
        "dct_done": False,
        # Per-step "Continue" gating flags. A step's Continue button is
        # enabled iff every prior step's flag here is True.
        "step_done": {
            "upload": False,
            "mapping": False,
            "groups": False,
            "outliers": False,
            # housekeeping completion is governed by hk_applied above
        },
    }


def _refresh_step_gates(state: dict, refs: dict) -> None:
    """Enable/disable each step's Continue button based on prior progress.

    Rule: a step's Continue is enabled iff every prior step has been
    completed. The very first step's Continue is enabled the moment its
    own work is done (e.g. file upload finishes).
    """
    sd = state["step_done"]
    ordered = ["upload", "mapping", "groups", "outliers"]

    def _set(name: str, enabled: bool) -> None:
        btn = refs.get(f"{name}_next_btn")
        if btn is None:
            return
        if enabled:
            btn.enable()
        else:
            btn.disable()

    # Upload: enabled once a file has been read.
    _set("upload", sd["upload"])
    # Each subsequent step's Continue is gated on every earlier flag.
    cumulative = sd["upload"]
    _set("mapping", cumulative)
    cumulative = cumulative and sd["mapping"]
    _set("groups", cumulative)
    cumulative = cumulative and sd["groups"]
    _set("outliers", cumulative)
    cumulative = cumulative and sd["outliers"]
    # Housekeeping Continue also requires hk_applied; handled in its builder.
    if "hk_next_btn" in refs:
        if cumulative and state.get("hk_applied"):
            refs["hk_next_btn"].enable()
        else:
            refs["hk_next_btn"].disable()
    # ΔCt run button: needs all five prior steps done.
    if "dct_run_btn" in refs:
        if cumulative and state.get("hk_applied"):
            refs["dct_run_btn"].enable()
        else:
            refs["dct_run_btn"].disable()


@ui.page("/")
def index() -> None:
    """Render the two-pane workflow."""
    state = _new_state()
    refs: dict[str, Any] = {}

    ui.query("body").classes("bg-slate-50")

    with ui.header(elevated=False).classes(
        "bg-white text-slate-900 items-center border-b border-slate-200 px-6"
    ):
        with ui.row().classes("items-center gap-3"):
            ui.icon("biotech", size="28px").classes("text-blue-600")
            ui.label("qPCR Analyzer").classes("text-xl font-semibold")
            ui.badge(f"v{__version__}", color="blue").classes(
                "text-white text-xs"
            )
            ui.label("· ΔCt & batch-aware ΔΔCt quantification").classes(
                "text-sm text-slate-500"
            )
        ui.space()

    with ui.row().classes("w-full no-wrap items-stretch gap-0 p-4"):
        # ── LEFT PANE: workflow stepper ───────────────────────────────────────
        with ui.column().classes(
            "w-2/5 min-w-[420px] max-w-[640px] bg-white border border-slate-200 "
            "rounded-lg shadow-sm p-4 mr-4"
        ):
            with ui.expansion("How to use this app", icon="help_outline").classes(
                "w-full"
            ):
                ui.markdown(
                    "**Workflow.** Work through the steps on the left. The right "
                    "panel updates in real time as you upload data, choose "
                    "settings, and run analyses.\n\n"
                    "1. **Upload** your qPCR file.\n"
                    "2. **Confirm column mapping** — the app auto-detects "
                    "Well / Target / Sample / Cq, plus optional Group and Batch "
                    "columns if your file already annotates them.\n"
                    "3. **Review groups & batches.** If your file did not "
                    "include them, assign them here. A sanity check verifies "
                    "each sample belongs to one group and one batch.\n"
                    "4. **Review outliers.** Tighten or relax the replicate "
                    "tolerance and confirm which wells should be excluded.\n"
                    "5. **Pick housekeeping gene(s).** Samples missing a valid "
                    "housekeeping Cq are surfaced for confirmation before "
                    "being skipped for that gene only. Click **Apply** to "
                    "lock in your selection.\n"
                    "6. **Run ΔCt.** Always required.\n"
                    "7. **Run ΔΔCt** (optional). Pick a reference group and "
                    "press the button — only enabled after ΔCt has run.\n\n"
                    "**ΔCt** = mean Cq(target) − mean Cq(housekeeping).  \n"
                    "**ΔΔCt** = ΔCt(sample) − mean ΔCt(reference group, batch).  \n"
                    "Relative expression = 2^(−ΔΔCt)."
                ).classes("text-sm text-slate-700")

            with ui.stepper().props("vertical flat header-nav").classes(
                "w-full"
            ) as stepper:
                refs["stepper"] = stepper
                _build_step_upload(state, refs, stepper)
                _build_step_mapping(state, refs, stepper)
                _build_step_groups(state, refs, stepper)
                _build_step_outliers(state, refs, stepper)
                _build_step_housekeeping(state, refs, stepper)
                _build_step_dct(state, refs, stepper)
                _build_step_ddct(state, refs, stepper)

        # ── RIGHT PANE: persistent output tabs ────────────────────────────────
        with ui.column().classes(
            "flex-grow bg-white border border-slate-200 rounded-lg shadow-sm p-4"
        ):
            ui.label("Output").classes("text-base font-semibold text-slate-800 mb-1")
            with ui.tabs().classes("w-full") as out_tabs:
                tab_summary = ui.tab("Summary", icon="summarize")
                tab_preview = ui.tab("Data preview", icon="table_view")
                tab_excluded = ui.tab("Excluded blocks", icon="rule")
                tab_hk_setup = ui.tab(
                    "Housekeeping & exclusions", icon="science"
                )
                tab_dct = ui.tab("ΔCt results", icon="show_chart")
                tab_ddct = ui.tab("ΔΔCt results", icon="bar_chart")
                tab_downloads = ui.tab("Downloads", icon="download")
            refs["out_tabs"] = out_tabs

            with ui.tab_panels(out_tabs, value=tab_summary).classes("w-full"):
                with ui.tab_panel(tab_summary):
                    refs["summary_panel"] = ui.column().classes("w-full gap-2")
                    with refs["summary_panel"]:
                        ui.label("Upload a file to see the dataset summary.").classes(
                            "text-sm text-slate-500"
                        )

                with ui.tab_panel(tab_preview):
                    refs["preview_panel"] = ui.column().classes("w-full gap-2")
                    with refs["preview_panel"]:
                        ui.label("Data preview will appear after column mapping.").classes(
                            "text-sm text-slate-500"
                        )

                with ui.tab_panel(tab_excluded):
                    refs["excluded_panel"] = ui.column().classes("w-full gap-2")
                    with refs["excluded_panel"]:
                        ui.label(
                            "Outlier review and exclusion blocks will appear here."
                        ).classes("text-sm text-slate-500")

                with ui.tab_panel(tab_hk_setup):
                    refs["hk_setup_panel"] = ui.column().classes(
                        "w-full gap-3"
                    )
                    with refs["hk_setup_panel"]:
                        ui.label(
                            "Housekeeping gene selection and the associated "
                            "sample / well exclusions will appear here once "
                            "you reach step 5."
                        ).classes("text-sm text-slate-500")

                with ui.tab_panel(tab_dct):
                    refs["dct_panel"] = ui.column().classes("w-full gap-3")
                    with refs["dct_panel"]:
                        ui.label("Run ΔCt to populate this tab.").classes(
                            "text-sm text-slate-500"
                        )

                with ui.tab_panel(tab_ddct):
                    refs["ddct_panel"] = ui.column().classes("w-full gap-3")
                    with refs["ddct_panel"]:
                        ui.label("Run ΔΔCt to populate this tab.").classes(
                            "text-sm text-slate-500"
                        )

                with ui.tab_panel(tab_downloads):
                    refs["downloads_panel"] = ui.column().classes("w-full gap-2")
                    with refs["downloads_panel"]:
                        ui.label(
                            "Downloads (Excel workbook + per-figure PNG) appear "
                            "after analyses run."
                        ).classes("text-sm text-slate-500")

    with ui.footer().classes(
        "bg-white border-t border-slate-200 text-slate-600 text-xs px-6 py-3 "
        "items-center justify-center"
    ):
        with ui.row().classes("items-center gap-2 flex-wrap justify-center"):
            ui.label(f"qPCR Analyzer v{__version__}").classes("font-medium")
            ui.label("·").classes("text-slate-400")
            ui.label("Created by Jielin Yang")
            ui.label("·").classes("text-slate-400")
            ui.label("MIT License")
            ui.label("·").classes("text-slate-400")
            with ui.link(
                target="https://github.com/j-y26/py_qpcr_analyzer",
                new_tab=True,
            ).classes(
                "flex items-center justify-center w-7 h-7 rounded-full "
                "text-slate-500 hover:text-white hover:bg-slate-800 "
                "transition-colors no-underline"
            ).tooltip("View source on GitHub"):
                ui.html(
                    '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" '
                    'viewBox="0 0 24 24" fill="currentColor" aria-label="GitHub">'
                    '<path d="M12 .5C5.73.5.67 5.56.67 11.83c0 5.02 3.25 9.27 7.76 '
                    '10.77.57.1.78-.25.78-.55 0-.27-.01-.99-.02-1.94-3.16.69-3.83-'
                    '1.52-3.83-1.52-.52-1.32-1.27-1.67-1.27-1.67-1.04-.71.08-.7.08-.7 '
                    '1.15.08 1.76 1.18 1.76 1.18 1.02 1.75 2.69 1.25 3.35.96.1-.74.4-'
                    '1.25.72-1.54-2.52-.29-5.18-1.26-5.18-5.62 0-1.24.45-2.26 1.18-'
                    '3.06-.12-.29-.51-1.45.11-3.02 0 0 .96-.31 3.15 1.17a10.9 10.9 0 0 '
                    '1 5.74 0c2.18-1.48 3.14-1.17 3.14-1.17.63 1.57.23 2.73.11 3.02.74.'
                    '8 1.18 1.82 1.18 3.06 0 4.37-2.67 5.32-5.21 5.61.41.36.78 1.06.78 '
                    '2.13 0 1.54-.01 2.78-.01 3.16 0 .31.21.66.79.55 4.51-1.5 7.75-5.75 '
                    '7.75-10.77C23.33 5.56 18.27.5 12 .5z"/></svg>'
                )

    # Initial gating: nothing has been done yet, so disable everything.
    _refresh_step_gates(state, refs)


# ── Step builders ────────────────────────────────────────────────────────────

def _build_step_upload(state: dict, refs: dict, stepper) -> None:
    with ui.step("upload", title="1. Upload data", icon="upload_file"):
        ui.markdown(
            "Drop in a qPCR results file. Supported: **xlsx, xls, csv, tsv, "
            "txt** (one row per well)."
        ).classes("text-slate-700")

        async def on_upload(e: events.UploadEventArguments) -> None:
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
            _render_initial_summary(state, refs)
            _render_mapping(state, refs)
            state["step_done"]["upload"] = True
            # A new file invalidates everything downstream.
            for k in ("mapping", "groups", "outliers"):
                state["step_done"][k] = False
            state["hk_applied"] = False
            _refresh_step_gates(state, refs)

        ui.upload(on_upload=on_upload, auto_upload=True, max_files=1).props(
            "accept=.xlsx,.xls,.csv,.tsv,.txt color=primary flat bordered"
        ).classes("w-full")
        refs["file_info"] = ui.label("").classes("text-sm text-slate-600")

        def _go_next() -> None:
            stepper.next()

        with ui.stepper_navigation():
            refs["upload_next_btn"] = (
                ui.button("Continue", on_click=_go_next)
                .props("color=primary unelevated")
            )


def _build_step_mapping(state: dict, refs: dict, stepper) -> None:
    with ui.step("mapping", title="2. Column mapping", icon="view_column"):
        ui.markdown(
            "Confirm the auto-detected column for each role. **Group** and "
            "**Batch** are optional — if your file does not annotate them, "
            "leave them as `(none)` and assign in the next step."
        ).classes("text-slate-700")
        refs["mapping_container"] = ui.column().classes("w-full gap-2")

        def go_to_groups() -> None:
            mapping = state["mapping"]
            if mapping is None:
                ui.notify("Upload a file first.", type="negative")
                return
            errs = mapping.validate()
            if errs:
                for msg in errs:
                    ui.notify(msg, type="negative")
                return
            std = apply_mapping(state["raw_df"], mapping)
            if std["Cq"].notna().sum() == 0:
                ui.notify(
                    "Selected Cq column has no numeric values.", type="negative"
                )
                return

            sg_errs = validate_sample_groups(std)
            if sg_errs:
                for msg in sg_errs:
                    ui.notify(msg, type="negative")
                return
            sb_errs = validate_sample_batches(std)
            if sb_errs:
                for msg in sb_errs:
                    ui.notify(msg, type="negative")
                return

            state["standardized"] = std
            state["_outliers_initialized"] = False
            state["excluded_wells"] = set()
            state["sample_excludes_per_hk"] = {}
            # Seed sample_batches from the file if a Batch column is present,
            # and pre-tick the multi-batch checkbox iff the file actually
            # spans more than one batch.
            if "Batch" in std.columns:
                state["sample_batches"] = (
                    std.groupby("Sample")["Batch"].first().astype(str).to_dict()
                )
                state["has_batches"] = (
                    len(set(state["sample_batches"].values())) > 1
                )
            else:
                state["sample_batches"] = {}
                state["has_batches"] = False
            state["reference_group"] = None
            state["dct_done"] = False
            state["dct_results"] = None
            state["ddct_results"] = None
            state["summary"] = summarize_dataset(std, filename=state["filename"])
            _render_full_summary(state, refs)
            _render_data_preview(state, refs)
            _render_groups(state, refs)
            state["step_done"]["mapping"] = True
            for k in ("groups", "outliers"):
                state["step_done"][k] = False
            state["hk_applied"] = False
            _refresh_step_gates(state, refs)
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            refs["mapping_next_btn"] = ui.button(
                "Continue", on_click=go_to_groups
            ).props("color=primary unelevated")


def _build_step_groups(state: dict, refs: dict, stepper) -> None:
    with ui.step("groups", title="3. Groups & batches", icon="groups"):
        ui.markdown(
            "Verify which **group** and **batch** each sample belongs to. "
            "If your file already includes the Group/Batch columns, the table "
            "is pre-filled and you can edit values inline. The ΔΔCt method "
            "normalises within each batch before merging — leave all batches "
            "as `batch_1` for a single experimental run."
        ).classes("text-slate-700")
        refs["groups_container"] = ui.column().classes("w-full gap-3")

        def go_to_outliers() -> None:
            if state.get("standardized") is None:
                ui.notify("Confirm column mapping first.", type="negative")
                return
            multi_batch = bool(state.get("has_batches"))
            # Collect inline edits from the table
            edited = refs.get("groups_rows") or []
            if edited:
                rows_df = pd.DataFrame(edited)
                std = state["standardized"].copy()
                # Apply group overrides; batch overrides only apply when the
                # multi-batch checkbox is on.
                group_map = dict(zip(rows_df["Sample"], rows_df["Group"]))
                std["Group"] = std["Sample"].map(group_map).fillna(std["Group"]).astype(str)
                if multi_batch:
                    batch_map = dict(zip(rows_df["Sample"], rows_df["Batch"]))
                    std["Batch"] = (
                        std["Sample"].map(batch_map).fillna("batch_1").astype(str)
                    )
                    state["sample_batches"] = {
                        str(s): str(b) for s, b in batch_map.items()
                    }
                else:
                    # Drop any Batch column carried in from the file and
                    # collapse the in-memory batch map to a single batch.
                    if "Batch" in std.columns:
                        std = std.drop(columns=["Batch"])
                    state["sample_batches"] = {}
                state["standardized"] = std
                state["summary"] = summarize_dataset(std, filename=state["filename"])
                _render_full_summary(state, refs)
                _render_data_preview(state, refs)

            sg_errs = validate_sample_groups(state["standardized"])
            for msg in sg_errs:
                ui.notify(msg, type="negative")
            if sg_errs:
                return
            if multi_batch:
                sb_errs = validate_sample_batches(state["standardized"])
                for msg in sb_errs:
                    ui.notify(msg, type="negative")
                if sb_errs:
                    return

            _render_outliers(state, refs)
            state["step_done"]["groups"] = True
            state["step_done"]["outliers"] = False
            state["hk_applied"] = False
            _refresh_step_gates(state, refs)
            # Show the user the right-pane *Excluded blocks* tab so they can
            # see flagged replicate blocks while picking outliers.
            refs["out_tabs"].set_value("Excluded blocks")
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            refs["groups_next_btn"] = ui.button(
                "Continue", on_click=go_to_outliers
            ).props("color=primary unelevated")


def _build_step_outliers(state: dict, refs: dict, stepper) -> None:
    with ui.step(
        "outliers",
        title="4. Outliers",
        icon="rule",
    ):
        with ui.expansion("Outlier rule", icon="info").classes("w-full"):
            ui.markdown(
                "Within each (Sample × Target) replicate set, the app keeps "
                "the tightest contiguous run whose Cq range is within the "
                "tolerance and flags the rest. NaN Cq values are always "
                "flagged. Excluded blocks float to the top of the right-pane "
                "*Excluded blocks* tab — replicates of any flagged or "
                "manually excluded well are shown together so you can see "
                "context, not just the offending well."
            ).classes("text-sm text-slate-600")
        refs["outlier_container"] = ui.column().classes("w-full gap-3")

        def go_to_housekeeping() -> None:
            if state.get("standardized") is None:
                ui.notify("Confirm column mapping first.", type="negative")
                return
            _render_housekeeping(state, refs)
            state["step_done"]["outliers"] = True
            state["hk_applied"] = False
            _refresh_step_gates(state, refs)
            # Show the user the right-pane housekeeping/exclusions tab so
            # they can see the live summary while configuring step 5.
            refs["out_tabs"].set_value("Housekeeping & exclusions")
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            refs["outliers_next_btn"] = ui.button(
                "Continue", on_click=go_to_housekeeping
            ).props("color=primary unelevated")


def _build_step_housekeeping(state: dict, refs: dict, stepper) -> None:
    with ui.step(
        "housekeeping",
        title="5. Housekeeping gene(s)",
        icon="science",
    ):
        ui.markdown(
            "Pick the housekeeping gene(s) to normalise against. ΔCt is "
            "computed *per housekeeping gene*, so picking more than one "
            "yields parallel sets of results. Samples missing a valid "
            "housekeeping Cq are surfaced for confirmation before being "
            "skipped — for that gene only.\n\n"
            "**Click Apply to lock in your selection** before continuing to "
            "the ΔCt step."
        ).classes("text-slate-700")

        refs["hk_container"] = ui.column().classes("w-full gap-3")

        def go_to_dct() -> None:
            if not state.get("ref_genes"):
                ui.notify("Pick at least one housekeeping gene.", type="negative")
                return
            if not state.get("hk_applied"):
                ui.notify(
                    "Click Apply to lock in your housekeeping gene selection.",
                    type="negative",
                )
                return
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            hk_next = ui.button("Continue", on_click=go_to_dct).props(
                "color=primary unelevated"
            )
            hk_next.disable()
            refs["hk_next_btn"] = hk_next


def _build_step_dct(state: dict, refs: dict, stepper) -> None:
    with ui.step("dct", title="6. Run ΔCt", icon="play_arrow"):
        ui.markdown(
            "ΔCt normalises each sample to its housekeeping gene Cq. "
            "No reference group is required at this step — the result is "
            "ready as soon as you click below.\n\n"
            "Samples flagged in the previous step as *missing housekeeping "
            "Cq* will be skipped for that gene only."
        ).classes("text-slate-700")
        refs["dct_status"] = ui.label("").classes("text-sm text-slate-600")

        # Loading dialog used to keep the user informed during compute.
        with ui.dialog().props("persistent") as dct_dialog, ui.card().classes(
            "items-center gap-2 p-6"
        ):
            ui.spinner(size="lg", color="primary")
            ui.label("Running ΔCt analysis…").classes(
                "text-base font-semibold text-slate-800"
            )
            ui.label(
                "Crunching mean Cq values, applying exclusions, and building "
                "result tables. This usually takes a few seconds."
            ).classes("text-xs text-slate-500 text-center max-w-xs")
        refs["dct_loading_dialog"] = dct_dialog

        async def run_dct() -> None:
            if state["standardized"] is None:
                ui.notify("Load and configure data first.", type="negative")
                return
            if not state.get("ref_genes"):
                ui.notify("Pick at least one housekeeping gene.", type="negative")
                return
            std = state["standardized"].copy()
            std["Excluded"] = std["Well"].isin(state["excluded_wells"])
            dct_dialog.open()
            try:
                # Yield to the event loop so the dialog actually paints
                # before we kick off the (synchronous) compute.
                await asyncio.sleep(0.05)
                dct = await asyncio.to_thread(
                    compute_delta_ct,
                    std,
                    list(state["ref_genes"]),
                    state["sample_excludes_per_hk"],
                )
            except Exception as ex:  # noqa: BLE001
                dct_dialog.close()
                ui.notify(f"ΔCt failed: {ex}", type="negative")
                return
            dct_dialog.close()
            state["dct_results"] = dct
            state["dct_done"] = True
            refs["dct_status"].set_text(
                f"ΔCt computed for {len(dct)} housekeeping gene(s)."
            )
            _render_dct_results(state, refs)
            _render_downloads(state, refs)
            _render_ddct_setup(state, refs)
            if "ddct_run_btn" in refs:
                refs["ddct_run_btn"].enable()
            refs["out_tabs"].set_value("ΔCt results")
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            refs["dct_run_btn"] = ui.button("Run ΔCt", on_click=run_dct).props(
                "color=primary unelevated icon=play_arrow"
            )
            refs["dct_run_btn"].disable()


def _build_step_ddct(state: dict, refs: dict, stepper) -> None:
    with ui.step("ddct", title="7. Run ΔΔCt (optional)", icon="bar_chart"):
        ui.markdown(
            "Pick a **reference group**. Within each batch, the mean ΔCt of "
            "the reference group's samples anchors the relative expression "
            "scale; ΔΔCt = ΔCt(sample) − ΔCt(reference, batch). This step is "
            "optional — ΔCt alone is enough for housekeeping-normalised "
            "comparisons."
        ).classes("text-slate-700")
        refs["ddct_container"] = ui.column().classes("w-full gap-2")
        with refs["ddct_container"]:
            ui.label(
                "Reference group, batch summary, and the run button appear "
                "here once ΔCt has succeeded."
            ).classes("text-xs text-slate-500")

        with ui.dialog().props("persistent") as ddct_dialog, ui.card().classes(
            "items-center gap-2 p-6"
        ):
            ui.spinner(size="lg", color="primary")
            ui.label("Running ΔΔCt analysis…").classes(
                "text-base font-semibold text-slate-800"
            )
            ui.label(
                "Anchoring the reference group within each batch and computing "
                "relative expression. This usually takes a few seconds."
            ).classes("text-xs text-slate-500 text-center max-w-xs")
        refs["ddct_loading_dialog"] = ddct_dialog

        async def run_ddct() -> None:
            if not state["dct_done"]:
                ui.notify("Run ΔCt first.", type="negative")
                return
            ref_sel = refs.get("ref_group_sel")
            ref_group = ref_sel.value if ref_sel else None
            if not ref_group:
                ui.notify("Pick a reference group.", type="negative")
                return
            state["reference_group"] = ref_group
            std = state["standardized"].copy()
            std["Excluded"] = std["Well"].isin(state["excluded_wells"])
            ddct_dialog.open()
            try:
                await asyncio.sleep(0.05)
                ddct = await asyncio.to_thread(
                    compute_delta_delta_ct,
                    std,
                    list(state["ref_genes"]),
                    ref_group,
                    state["sample_batches"] or None,
                    state["sample_excludes_per_hk"],
                )
            except Exception as ex:  # noqa: BLE001
                ddct_dialog.close()
                ui.notify(f"ΔΔCt failed: {ex}", type="negative")
                return
            ddct_dialog.close()
            state["ddct_results"] = ddct
            _render_ddct_results(state, refs)
            _render_downloads(state, refs)
            refs["out_tabs"].set_value("ΔΔCt results")

        refs["ddct_run_fn"] = run_ddct

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            run_btn = ui.button("Run ΔΔCt", on_click=run_ddct).props(
                "color=primary unelevated icon=play_arrow"
            )
            run_btn.disable()
            refs["ddct_run_btn"] = run_btn


def _render_ddct_setup(state: dict, refs: dict) -> None:
    """Populate step 6 with the reference-group select + batch summary.

    Called from the ΔCt run handler so the reference-group dropdown is
    pre-filled with the actual groups present in the standardised data.
    """
    container = refs.get("ddct_container")
    if container is None:
        return
    container.clear()
    df = state["standardized"]
    groups = list(pd.unique(df["Group"].astype(str))) if df is not None else []
    batches = sorted(
        set((state.get("sample_batches") or {}).values()) or {"batch_1"}
    )
    with container:
        ui.label("ΔΔCt Configuration").classes("text-base font-semibold")
        with ui.row().classes("w-full items-end gap-3 flex-wrap"):
            sel = (
                ui.select(
                    options=groups,
                    label="Reference group (ΔΔCt anchor)",
                    value=state.get("reference_group") or (groups[0] if groups else None),
                )
                .classes("min-w-64")
                .props("outlined dense")
            )
            refs["ref_group_sel"] = sel
        ui.label(
            f"{len(batches)} batch(es): {', '.join(batches)}"
        ).classes("text-xs text-slate-500")
        ui.label(
            "Each batch must contain at least one sample from the reference "
            "group for every measured target. Click ‘Run ΔΔCt’ below when "
            "ready."
        ).classes("text-xs text-slate-500")


# ── Render helpers ────────────────────────────────────────────────────────────

def _render_initial_summary(state: dict, refs: dict) -> None:
    """Right-pane summary using the raw uploaded DataFrame.

    A best-effort summary that renders before the user confirms column
    mapping. Once mapping is applied a richer summary replaces this one.
    """
    panel = refs["summary_panel"]
    panel.clear()
    df = state["raw_df"]
    n_rows, n_cols = df.shape
    with panel:
        ui.label("Initial dataset preview").classes("text-base font-semibold")
        ui.label(
            f"Loaded {state['filename']} · {n_rows} rows × {n_cols} columns"
        ).classes("text-sm text-slate-700")
        ui.label(
            "Confirm column mapping in step 2 to see the full sample / "
            "target / group breakdown."
        ).classes("text-xs text-slate-500")


def _render_full_summary(state: dict, refs: dict) -> None:
    panel = refs["summary_panel"]
    panel.clear()
    s = state["summary"]
    df = state.get("standardized")
    show_batches = bool(state.get("has_batches"))
    with panel:
        ui.label("Dataset summary").classes(
            "text-lg font-semibold text-slate-800"
        )
        with ui.row().classes("w-full gap-6 flex-wrap text-sm text-slate-700"):
            ui.label(f"Analysis: {s['analysis_time']}")
            ui.label(f"File: {s['filename']}")
        with ui.row().classes("w-full gap-2 flex-wrap"):
            _stat_chip("Wells", s["n_wells"])
            _stat_chip("NA Cq wells", s["n_na_wells"], warn=s["n_na_wells"] > 0)
            _stat_chip("Replicate blocks", s["n_replicate_blocks"])
            _stat_chip(
                "Replicates / block",
                f"{s['replicates_min']}–{s['replicates_max']}",
            )

        # ── Samples grouped by Group ────────────────────────────────────
        if df is not None and "Group" in df.columns:
            samples_by_group: dict[str, list[str]] = {}
            for g in group_order(df):
                samples_by_group[g] = list(
                    pd.unique(df.loc[df["Group"] == g, "Sample"].astype(str))
                )
            with ui.expansion(
                f"Samples: {s['n_samples']}",
                icon="people_alt",
            ).classes(
                "w-full rounded-md border border-slate-200 bg-white shadow-sm"
            ):
                with ui.column().classes("w-full gap-3 p-1"):
                    for i, (g, samples) in enumerate(samples_by_group.items()):
                        accent = _ACCENTS[i % len(_ACCENTS)]
                        with ui.column().classes("w-full gap-1"):
                            with ui.row().classes("items-center gap-2"):
                                ui.label(g).classes(
                                    f"text-sm font-semibold {accent['text']}"
                                )
                                ui.badge(
                                    str(len(samples)),
                                    color=accent["badge"],
                                ).classes("text-white")
                            with ui.row().classes("w-full gap-1 flex-wrap"):
                                for v in samples:
                                    _pill(v, accent)
        else:
            _summary_section_pills("Samples", s["n_samples"], s["samples"], 0)

        _summary_section_pills(
            "Targets", s["n_targets"], s["targets"], accent_idx=2
        )
        _summary_section_pills(
            "Groups", s["n_groups"], s["groups"], accent_idx=4
        )
        if show_batches and s["has_batch_column"]:
            _summary_section_pills(
                "Batches", s["n_batches"], s["batches"], accent_idx=5
            )


# A small palette of soft tints used by the summary pills/section headings.
_ACCENTS: list[dict[str, str]] = [
    {"text": "text-blue-700",   "bg": "bg-blue-50",    "border": "border-blue-200",    "badge": "blue"},
    {"text": "text-emerald-700","bg": "bg-emerald-50", "border": "border-emerald-200", "badge": "green"},
    {"text": "text-violet-700", "bg": "bg-violet-50",  "border": "border-violet-200",  "badge": "purple"},
    {"text": "text-amber-700",  "bg": "bg-amber-50",   "border": "border-amber-200",   "badge": "amber"},
    {"text": "text-rose-700",   "bg": "bg-rose-50",    "border": "border-rose-200",    "badge": "red"},
    {"text": "text-cyan-700",   "bg": "bg-cyan-50",    "border": "border-cyan-200",    "badge": "cyan"},
]


def _pill(text: str, accent: dict[str, str]) -> None:
    """Render a single soft-tinted chip used in dataset-summary sections."""
    with ui.element("div").classes(
        f"px-2.5 py-1 rounded-full text-xs {accent['bg']} "
        f"border {accent['border']} {accent['text']} font-medium"
    ):
        ui.label(str(text))


def _stat_chip(label: str, value, warn: bool = False) -> None:
    color = "amber" if warn else "blue"
    value_str = "" if value is None else str(value)
    with ui.element("div").classes(
        f"px-3 py-1 rounded-full bg-{color}-50 border border-{color}-200 "
        "text-sm flex items-center gap-1"
    ):
        if value_str:
            ui.label(value_str).classes(f"font-semibold text-{color}-700")
        ui.label(label).classes(f"text-{color}-700 font-medium")


def _summary_section_pills(
    title: str,
    count: int,
    values: list[str],
    accent_idx: int,
) -> None:
    accent = _ACCENTS[accent_idx % len(_ACCENTS)]
    with ui.expansion(
        f"{title}: {count}",
        icon="list",
    ).classes(
        "w-full rounded-md border border-slate-200 bg-white shadow-sm"
    ):
        if not values:
            ui.label("(none)").classes("text-sm text-slate-500")
            return
        with ui.row().classes("w-full gap-1 flex-wrap p-1"):
            for v in values:
                _pill(v, accent)


def _summary_section(title: str, count: int, values: list[str]) -> None:
    """Backwards-compatible plain-grey list (kept for any future reuse)."""
    with ui.expansion(f"{title}: {count}", icon="list").classes(
        "w-full bg-slate-50 rounded"
    ):
        if not values:
            ui.label("(none)").classes("text-sm text-slate-500")
            return
        with ui.row().classes("w-full gap-1 flex-wrap"):
            for v in values:
                ui.badge(str(v), color="grey").classes("text-white")


def _render_data_preview(state: dict, refs: dict) -> None:
    panel = refs["preview_panel"]
    panel.clear()
    df = state["standardized"]
    if df is None:
        return
    sorted_df = sort_wells(df)
    show_batches = bool(state.get("has_batches"))
    candidate_cols = ["Well", "Target", "Group", "Sample", "Batch", "Cq"]
    cols = [
        c for c in candidate_cols
        if c in sorted_df.columns and (c != "Batch" or show_batches)
    ]
    rows = [
        {c: (None if pd.isna(r[c]) else r[c]) for c in cols}
        for _, r in sorted_df[cols].iterrows()
    ]
    columns = [
        {"name": c, "label": c, "field": c, "align": "left", "sortable": True}
        for c in cols
    ]
    with panel:
        ui.label(
            "Standardised data — sorted by Target (file order) → Group → "
            "Sample → Well."
        ).classes("text-sm text-slate-600")
        ui.table(columns=columns, rows=rows, row_key="Well", pagination=20).classes(
            "w-full"
        )


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
        df[["Sample", "Group"]].drop_duplicates()
        .sort_values(["Group", "Sample"]).reset_index(drop=True)
    )

    has_batch_col = "Batch" in df.columns
    if has_batch_col:
        sample_batch = df.groupby("Sample")["Batch"].first().astype(str).to_dict()
    else:
        sample_batch = {}

    rows = [
        {
            "Sample": r.Sample,
            "Group": r.Group,
            "Batch": sample_batch.get(r.Sample, "batch_1"),
        }
        for r in samples_meta.itertuples()
    ]
    refs["groups_rows"] = rows
    initial_groups = sorted({str(r["Group"]) for r in rows})
    initial_batches = sorted({str(r["Batch"]) for r in rows})

    with container:
        ui.label(
            f"{len(samples_meta)} unique sample(s). "
            + (
                "Batch column detected — pre-filled from file."
                if has_batch_col
                else "No Batch column in file."
            )
        ).classes("text-sm text-slate-600")

        # ── Multi-batch toggle ────────────────────────────────────────────
        batch_toggle = ui.checkbox(
            "Samples are from different batches (enable Batch column)",
            value=bool(state.get("has_batches")),
        ).classes("text-sm")
        ui.label(
            "When unchecked, every sample is treated as one batch and the "
            "Batch column is excluded from the data preview, summary, and "
            "exported files."
        ).classes("text-xs text-slate-500")

        ui.label(
            "Group / Batch cells accept either a pick from the dropdown or a "
            "freshly-typed value (press Enter to commit). New values join the "
            "dropdown for subsequent rows."
        ).classes("text-xs text-slate-500")

        with ui.card().classes("w-full border border-slate-200 shadow-none p-3"):
            grid_classes = (
                "grid grid-cols-3 gap-x-3 gap-y-2 items-center w-full"
            )
            with ui.element("div").classes(
                grid_classes + " text-xs font-semibold text-slate-600"
            ):
                ui.label("Sample")
                ui.label("Group")
                ui.label("Batch")

            group_selects: list[ui.select] = []
            batch_selects: list[ui.select] = []

            def _refresh_group_options() -> None:
                opts = sorted(
                    {str(s.value) for s in group_selects if s.value}
                    | set(initial_groups)
                )
                for s in group_selects:
                    s.options = opts
                    s.update()

            def _refresh_batch_options() -> None:
                opts = sorted(
                    {str(s.value) for s in batch_selects if s.value}
                    | set(initial_batches)
                )
                for s in batch_selects:
                    s.options = opts
                    s.update()

            with ui.element("div").classes(grid_classes):
                for row in rows:
                    ui.label(row["Sample"]).classes("text-sm truncate")

                    g_sel = (
                        ui.select(
                            options=initial_groups,
                            value=row["Group"],
                            with_input=True,
                            new_value_mode="add-unique",
                        )
                        .classes("w-full")
                        .props("outlined dense")
                    )

                    def _on_group(_e, row=row, sel=None) -> None:
                        v = (str(sel.value).strip() if sel.value else "")
                        row["Group"] = v or "unassigned"
                        _refresh_group_options()

                    g_sel.on_value_change(
                        lambda e, row=row, sel=g_sel: _on_group(e, row=row, sel=sel)
                    )
                    group_selects.append(g_sel)

                    b_sel = (
                        ui.select(
                            options=initial_batches,
                            value=row["Batch"],
                            with_input=True,
                            new_value_mode="add-unique",
                        )
                        .classes("w-full")
                        .props("outlined dense")
                    )

                    def _on_batch(_e, row=row, sel=None) -> None:
                        v = (str(sel.value).strip() if sel.value else "")
                        row["Batch"] = v or "batch_1"
                        state["sample_batches"][str(row["Sample"])] = row["Batch"]
                        _refresh_batch_options()

                    b_sel.on_value_change(
                        lambda e, row=row, sel=b_sel: _on_batch(e, row=row, sel=sel)
                    )
                    batch_selects.append(b_sel)

        def _apply_batch_enabled(enabled: bool) -> None:
            for s in batch_selects:
                if enabled:
                    s.enable()
                else:
                    s.disable()

        def _on_batch_toggle(_e=None) -> None:
            state["has_batches"] = bool(batch_toggle.value)
            _apply_batch_enabled(state["has_batches"])

        batch_toggle.on_value_change(_on_batch_toggle)
        # Apply initial enabled state right after the dropdowns exist.
        _apply_batch_enabled(bool(state.get("has_batches")))

        ui.label(
            "Tip: leave the checkbox unchecked for a single experimental run."
        ).classes("text-xs text-slate-400")


def _render_outliers(state: dict, refs: dict) -> None:
    container = refs["outlier_container"]
    container.clear()

    with container:
        with ui.row().classes("w-full items-end gap-4 flex-wrap"):
            tol_input = ui.number(
                label="Replicate tolerance (cycles)",
                value=state["tolerance"],
                min=0.1,
                step=0.1,
                format="%.2f",
            ).classes("w-48").props("outlined dense")

            def _apply_tol() -> None:
                try:
                    new_tol = float(tol_input.value)
                except (TypeError, ValueError):
                    new_tol = 1.0
                if new_tol != state["tolerance"]:
                    state["tolerance"] = new_tol
                    state["_outliers_initialized"] = False
                _rebuild_outlier_view(state, refs)

            ui.button("Apply tolerance", on_click=_apply_tol).props(
                "color=primary outline"
            )

        refs["outlier_summary_slot"] = ui.label("").classes("text-sm text-slate-600")
        refs["excluded_select_slot"] = ui.column().classes("w-full")

        _rebuild_outlier_view(state, refs)


def _rebuild_outlier_view(state: dict, refs: dict) -> None:
    flagged = mark_outliers(state["standardized"], tolerance=state["tolerance"])
    state["flagged"] = flagged
    auto_wells = flagged.loc[flagged["Outlier"] & flagged["Cq"].notna(), "Well"].tolist()
    nan_wells = flagged.loc[flagged["Cq"].isna(), "Well"].tolist()
    if not state.get("_outliers_initialized"):
        state["excluded_wells"] = set(auto_wells) | set(nan_wells)
        state["_outliers_initialized"] = True

    sel_slot = refs["excluded_select_slot"]
    sel_slot.clear()
    targets = target_order(flagged)
    sorted_wells = sort_wells(flagged, targets=targets)["Well"].tolist()
    with sel_slot:
        excluded_select = (
            ui.select(
                options=sorted_wells,
                value=sorted(state["excluded_wells"], key=sorted_wells.index)
                if state["excluded_wells"]
                else [],
                multiple=True,
                label="Wells to exclude from analysis",
            )
            .classes("w-full")
            .props("outlined dense use-chips clearable")
        )

        def _on_change(_e) -> None:
            state["excluded_wells"] = set(excluded_select.value or [])
            refs["outlier_summary_slot"].set_text(_excluded_summary(state))
            _render_excluded_blocks(state, refs)
            # Excluded wells affect which samples lack a valid housekeeping
            # Cq, so the user must re-Apply on step 5 if they had already
            # confirmed.
            state["hk_applied"] = False
            _refresh_step_gates(state, refs)
            _render_excluded_samples_panel(state, refs)

        excluded_select.on_value_change(_on_change)

    refs["outlier_summary_slot"].set_text(_excluded_summary(state))
    _render_excluded_blocks(state, refs)
    _render_excluded_samples_panel(state, refs)


def _render_housekeeping(state: dict, refs: dict) -> None:
    """Render the housekeeping-gene selection step (gene picker + per-gene exclusion)."""
    container = refs.get("hk_container")
    if container is None:
        return
    container.clear()

    df = state["standardized"]
    targets = target_order(df)
    # Honour file-appearance order for sample names — biological labels
    # like donor_3 / donor_10 must not be lexicographically reordered.
    all_samples = sample_order(df)
    std = df.copy()
    std["Excluded"] = std["Well"].isin(state["excluded_wells"])

    excludes = state.setdefault("sample_excludes_per_hk", {})
    # Map of hk_gene → ui.select for the per-gene dropdowns. Mutated on
    # every render of the per-HK block so the global-exclude handler can
    # mirror its selection into each per-gene dropdown.
    per_hk_selects: dict[str, Any] = {}

    def _mark_dirty() -> None:
        """User changed something — they must click Apply again to continue."""
        state["hk_applied"] = False
        if "hk_next_btn" in refs:
            refs["hk_next_btn"].disable()
        if "hk_status_label" in refs:
            refs["hk_status_label"].set_text(
                "Selection changed — click Apply to lock it in."
            )
            refs["hk_status_label"].classes(
                replace="text-sm text-amber-700"
            )
        _render_excluded_samples_panel(state, refs)

    with container:
        ref_gene_sel = ui.select(
            options=targets,
            label="Housekeeping gene(s)",
            multiple=True,
            value=state["ref_genes"],
        ).classes("w-full").props("outlined dense use-chips")

        def _on_ref_change(_e) -> None:
            _mark_dirty()
            _render_per_hk_block()

        ref_gene_sel.on_value_change(_on_ref_change)

        per_hk_slot = ui.column().classes("w-full gap-2")

        def _render_per_hk_block() -> None:
            per_hk_slot.clear()
            per_hk_selects.clear()
            picked = list(ref_gene_sel.value or [])
            if not picked:
                with per_hk_slot:
                    ui.label(
                        "Pick at least one housekeeping gene to configure "
                        "per-gene sample exclusion."
                    ).classes("text-xs text-slate-500")
                return

            missing = samples_missing_hk(std, picked)
            global_excl = set(excludes.get("*", set()))
            with per_hk_slot:
                with ui.expansion(
                    "Why exclude samples per housekeeping gene?", icon="info"
                ).classes("w-full"):
                    ui.markdown(
                        "ΔCt is computed *per housekeeping gene*. If a sample "
                        "has no valid Cq for one gene, it can still be "
                        "analysed against another. Excluding a sample for a "
                        "single gene keeps it in the rest of the analysis "
                        "instead of breaking the whole run."
                    ).classes("text-sm text-slate-600")

                for hk in picked:
                    auto_excl = set(missing.get(hk, []))
                    # Pre-populate per-gene selection with auto-flagged
                    # samples *and* anything chosen in "Exclude entirely".
                    current = (
                        set(excludes.get(hk, set())) | auto_excl | global_excl
                    )
                    excludes[hk] = current

                    with ui.card().classes(
                        "w-full border border-slate-200 shadow-none"
                    ):
                        ui.label(f"Housekeeping: {hk}").classes("font-semibold")
                        if auto_excl:
                            ui.label(
                                "Auto-excluded (no valid housekeeping Cq): "
                                + ", ".join(sorted(auto_excl))
                            ).classes("text-xs text-amber-700")
                        sel = (
                            ui.select(
                                options=all_samples,
                                value=sorted(current),
                                multiple=True,
                                label=f"Samples to exclude for {hk}",
                            )
                            .classes("w-full")
                            .props("outlined dense use-chips clearable")
                        )
                        per_hk_selects[hk] = sel

                        def _on_excl_change(_e, hk=hk, sel=sel) -> None:
                            excludes[hk] = set(sel.value or [])
                            _mark_dirty()

                        sel.on_value_change(_on_excl_change)

                with ui.card().classes(
                    "w-full border border-slate-200 shadow-none"
                ):
                    ui.label(
                        "Exclude entirely (all housekeeping genes)"
                    ).classes("font-semibold")
                    ui.label(
                        "Samples picked here are also added to every "
                        "per-gene dropdown above."
                    ).classes("text-xs text-slate-500")
                    global_sel = (
                        ui.select(
                            options=all_samples,
                            value=sorted(excludes.get("*", set())),
                            multiple=True,
                            label="Samples to exclude from every analysis",
                        )
                        .classes("w-full")
                        .props("outlined dense use-chips clearable")
                    )

                    def _on_global(_e, sel=global_sel) -> None:
                        new_global = set(sel.value or [])
                        prev_global = set(excludes.get("*", set()))
                        added = new_global - prev_global
                        removed = prev_global - new_global
                        excludes["*"] = new_global
                        # Mirror into per-gene dropdowns: add new, drop those
                        # only present because they were globally excluded.
                        for hk, hk_sel in per_hk_selects.items():
                            cur = set(excludes.get(hk, set()))
                            cur |= added
                            cur -= removed
                            excludes[hk] = cur
                            hk_sel.set_value(sorted(cur))
                        _mark_dirty()

                    global_sel.on_value_change(_on_global)

        def _apply() -> None:
            picked = list(ref_gene_sel.value or [])
            if not picked:
                ui.notify(
                    "Pick at least one housekeeping gene before applying.",
                    type="negative",
                )
                return
            state["ref_genes"] = picked
            # Drop stale per-gene exclusions for genes the user removed.
            for stale in [k for k in excludes if k != "*" and k not in picked]:
                excludes.pop(stale, None)
            state["hk_applied"] = True
            label_txt = (
                f"Applied: {', '.join(picked)}. You can now continue to ΔCt."
            )
            refs["hk_status_label"].set_text(label_txt)
            refs["hk_status_label"].classes(replace="text-sm text-emerald-700")
            _refresh_step_gates(state, refs)
            _render_excluded_samples_panel(state, refs)
            ui.notify("Housekeeping gene selection applied.", type="positive")

        with ui.row().classes("items-center gap-3"):
            ui.button("Apply", on_click=_apply).props(
                "color=primary unelevated icon=check"
            )
            refs["hk_status_label"] = ui.label(
                "Pick at least one housekeeping gene and click Apply."
            ).classes("text-sm text-slate-600")

        _render_per_hk_block()
        _render_excluded_samples_panel(state, refs)


def _status_pill(applied: bool) -> None:
    """Compact 'applied' / 'pending' pill, used in the housekeeping panel."""
    if applied:
        cls = "bg-emerald-50 border-emerald-200 text-emerald-700"
        text = "applied"
    else:
        cls = "bg-amber-50 border-amber-200 text-amber-700"
        text = "pending Apply"
    with ui.element("div").classes(
        f"px-3 py-1 rounded-full border text-sm flex items-center gap-2 {cls}"
    ):
        ui.icon("check_circle" if applied else "schedule").classes(
            "text-base"
        )
        ui.label("Status").classes("font-medium")
        ui.label(text).classes("font-semibold")


def _render_excluded_samples_panel(state: dict, refs: dict) -> None:
    """Right-pane *Housekeeping & exclusions* tab.

    Three side-by-side sections:
      a) Housekeeping gene(s) currently selected.
      b) Sample exclusion summary — global + per-gene exclusions.
      c) Well exclusion summary — outlier-flagged / manually-excluded wells.
    """
    panel = refs.get("hk_setup_panel")
    if panel is None:
        return
    panel.clear()

    df = state.get("standardized")
    excludes = state.get("sample_excludes_per_hk", {})
    excluded_wells = state.get("excluded_wells", set())
    ref_genes = list(state.get("ref_genes") or [])

    order_index: dict[str, int] = {}
    if df is not None:
        order_index = {s: i for i, s in enumerate(sample_order(df))}

    def _by_file_order(samples) -> list[str]:
        return sorted(samples, key=lambda s: order_index.get(str(s), 1 << 30))

    well_exclusions: dict[str, set[str]] = {}
    if df is not None and excluded_wells:
        ex = df[df["Well"].astype(str).isin(excluded_wells)]
        for _, row in ex.iterrows():
            well_exclusions.setdefault(str(row["Sample"]), set()).add(
                str(row["Target"])
            )

    global_excl = _by_file_order(set(excludes.get("*", set())))

    with panel:
        if df is None:
            ui.label("No data loaded yet.").classes("text-sm text-slate-500")
            return

        ui.label(
            "Configure ΔCt analysis — review your housekeeping pick, "
            "sample-level exclusions, and well-level exclusions side by side."
        ).classes("text-sm text-slate-600")

        # Three parallel section cards laid out as a responsive grid.
        with ui.row().classes(
            "w-full no-wrap items-stretch gap-3"
        ):
            # ── (a) Housekeeping gene(s) ─────────────────────────────────
            with ui.card().classes(
                "flex-1 min-w-[220px] border border-blue-200 bg-blue-50 shadow-none"
            ):
                with ui.row().classes("items-center gap-2"):
                    ui.icon("science").classes("text-blue-700")
                    ui.label("Housekeeping gene(s)").classes(
                        "text-base font-semibold text-blue-900"
                    )
                _status_pill(bool(state.get("hk_applied")))
                if not ref_genes:
                    ui.label(
                        "No housekeeping gene picked yet."
                    ).classes("text-sm text-slate-600")
                else:
                    ui.label(f"{len(ref_genes)} selected").classes(
                        "text-xs text-blue-700"
                    )
                    with ui.row().classes("w-full gap-1 flex-wrap"):
                        for g in ref_genes:
                            _pill(g, _ACCENTS[0])

            # ── (b) Sample exclusion summary ────────────────────────────
            with ui.card().classes(
                "flex-1 min-w-[260px] border border-slate-200 shadow-none"
            ):
                with ui.row().classes("items-center gap-2"):
                    ui.icon("person_off").classes("text-rose-700")
                    ui.label("Sample exclusion summary").classes(
                        "text-base font-semibold text-slate-800"
                    )
                ui.label(
                    "Live view of every sample-level exclusion currently in "
                    "effect."
                ).classes("text-xs text-slate-500")
                with ui.row().classes("w-full gap-2 flex-wrap"):
                    _stat_chip(
                        "Excluded entirely",
                        len(global_excl),
                        warn=len(global_excl) > 0,
                    )
                    _stat_chip(
                        "Housekeeping gene(s) selected",
                        len(ref_genes) if ref_genes else 0,
                        warn=not ref_genes,
                    )

                if global_excl:
                    with ui.element("div").classes(
                        "w-full mt-1 p-2 rounded border border-rose-200 bg-rose-50"
                    ):
                        ui.label(
                            "Excluded entirely (all housekeeping genes)"
                        ).classes("font-semibold text-rose-800 text-sm")
                        ui.label(
                            f"{len(global_excl)} sample(s): "
                            + ", ".join(global_excl)
                        ).classes("text-sm text-rose-700")

                if not ref_genes:
                    ui.label(
                        "Pick housekeeping gene(s) in step 5 to see per-gene "
                        "exclusions here."
                    ).classes("text-sm text-slate-500")
                else:
                    for hk in ref_genes:
                        hk_excl = _by_file_order(set(excludes.get(hk, set())))
                        only_hk = [s for s in hk_excl if s not in global_excl]
                        with ui.element("div").classes(
                            "w-full mt-1 p-2 rounded border "
                            + (
                                "border-amber-200 bg-amber-50"
                                if hk_excl
                                else "border-slate-200 bg-white"
                            )
                        ):
                            ui.label(f"Housekeeping: {hk}").classes(
                                "font-semibold text-sm"
                            )
                            if not hk_excl:
                                ui.label(
                                    "No samples excluded for this gene."
                                ).classes("text-xs text-slate-500")
                            else:
                                ui.label(
                                    f"{len(hk_excl)} sample(s) excluded: "
                                    + ", ".join(hk_excl)
                                ).classes("text-sm text-amber-800")
                                if only_hk:
                                    ui.label(
                                        "Excluded only for this gene: "
                                        + ", ".join(only_hk)
                                    ).classes("text-xs text-amber-700")

            # ── (c) Well exclusion summary ──────────────────────────────
            with ui.card().classes(
                "flex-1 min-w-[260px] border border-slate-200 shadow-none"
            ):
                with ui.row().classes("items-center gap-2"):
                    ui.icon("rule").classes("text-amber-700")
                    ui.label("Well exclusion summary").classes(
                        "text-base font-semibold text-slate-800"
                    )
                ui.label(
                    "Individual replicate wells flagged or manually excluded "
                    "in step 4. These rows are dropped from mean-Cq, but the "
                    "samples themselves remain unless also excluded above."
                ).classes("text-xs text-slate-500")
                with ui.row().classes("w-full gap-2 flex-wrap"):
                    _stat_chip(
                        "Excluded wells",
                        len(excluded_wells),
                        warn=len(excluded_wells) > 0,
                    )
                    _stat_chip(
                        "Affected samples",
                        len(well_exclusions),
                        warn=len(well_exclusions) > 0,
                    )
                if not excluded_wells:
                    ui.label("No wells excluded.").classes(
                        "text-sm text-slate-500"
                    )
                else:
                    with ui.expansion(
                        f"Wells by sample × target ({len(excluded_wells)} well(s))",
                        icon="science",
                    ).classes("w-full"):
                        ui.label(
                            "Excluded wells listed by sample × target so you "
                            "can spot accidental drops."
                        ).classes("text-xs text-slate-500")
                        for sample in _by_file_order(well_exclusions.keys()):
                            targets_set = well_exclusions[sample]
                            ordered_targets = [
                                t for t in target_order(df) if t in targets_set
                            ]
                            ui.label(
                                f"{sample}: " + ", ".join(ordered_targets)
                            ).classes("text-xs text-slate-700")


def _render_excluded_blocks(state: dict, refs: dict) -> None:
    """Right-pane *Excluded blocks* tab — sample × target replicate groups."""
    panel = refs["excluded_panel"]
    panel.clear()
    flagged = state.get("flagged")
    if flagged is None:
        return
    blocks = build_blocks(
        flagged, state["excluded_wells"], targets=target_order(flagged)
    )
    bad = [b for b in blocks if b["has_exclusion"]]
    good = [b for b in blocks if not b["has_exclusion"]]

    with panel:
        ui.label(
            "Replicate blocks containing excluded or flagged wells appear "
            "first. Each block lists every replicate so you can judge context, "
            "not just the offending well."
        ).classes("text-sm text-slate-600")
        ui.label(_excluded_summary(state)).classes("text-sm text-slate-700")

        with ui.row().classes("w-full gap-2"):
            _stat_chip(
                "Blocks with exclusion",
                len(bad),
                warn=len(bad) > 0,
            )
            _stat_chip("Clean blocks", len(good))

        if bad:
            ui.label("Exclusion-bearing blocks").classes(
                "text-base font-semibold mt-2"
            )
            for blk in bad:
                _render_block_card(blk)

        if good:
            ui.label(f"{len(good)} clean replicate block(s)").classes(
                "text-base font-semibold mt-2"
            )
            ui.label(
                "Each clean block is collapsed by default — click to inspect "
                "individual replicate Cq values."
            ).classes("text-xs text-slate-500")
            for blk in good:
                _render_block_expansion(blk)


def _render_block_card(blk: dict) -> None:
    """Card shown for blocks containing exclusions — always-open detail."""
    with ui.card().classes(
        "w-full border border-amber-300 bg-amber-50 shadow-none"
    ):
        _render_block_header(blk)
        _render_block_replicate_table(blk)


def _render_block_expansion(blk: dict) -> None:
    """Individually-expandable, collapsed-by-default card for clean blocks."""
    title = f"{blk['target']} · {blk['sample']} ({blk['group']})"
    caption = f"{len(blk['replicates'])} replicate(s) ok"
    with ui.expansion(title, caption=caption, icon="check").classes(
        "w-full bg-slate-50 rounded border border-slate-200"
    ):
        _render_block_replicate_table(blk)


def _render_block_header(blk: dict) -> None:
    with ui.row().classes("items-center gap-3 w-full"):
        ui.label(blk["target"]).classes("font-semibold")
        ui.label("·").classes("text-slate-400")
        ui.label(f"{blk['sample']} ({blk['group']})").classes("text-sm")
        ui.space()
        if blk["n_excluded"] or blk["n_flagged"]:
            ui.label(
                f"{blk['n_excluded']} excluded · {blk['n_flagged']} flagged "
                f"of {len(blk['replicates'])}"
            ).classes("text-xs text-amber-700")
        else:
            ui.label(f"{len(blk['replicates'])} replicate(s) ok").classes(
                "text-xs text-slate-500"
            )


def _render_block_replicate_table(blk: dict) -> None:
    cols = [
        {"name": "well", "label": "Well", "field": "well", "align": "left"},
        {"name": "cq", "label": "Cq", "field": "cq", "align": "left"},
        {"name": "status", "label": "Status", "field": "status", "align": "left"},
    ]
    rows = [
        {
            "well": r["well"],
            "cq": "—" if r["cq"] is None else f"{r['cq']:.3f}",
            "status": (
                "excluded" if r["excluded"]
                else ("NA" if r["is_nan"] else ("outlier" if r["outlier"] else "ok"))
            ),
        }
        for r in blk["replicates"]
    ]
    ui.table(columns=cols, rows=rows, row_key="well").classes("w-full")


def _excluded_summary(state: dict) -> str:
    n_total = 0 if state["flagged"] is None else len(state["flagged"])
    n_excl = len(state["excluded_wells"])
    return f"{n_excl} of {n_total} well(s) marked for exclusion."


def _render_dct_results(state: dict, refs: dict) -> None:
    panel = refs["dct_panel"]
    panel.clear()
    results = state["dct_results"] or {}
    color_map = _color_map_for_results(state["standardized"], results)
    with panel:
        ui.markdown(
            "**Relative expression vs housekeeping** = 2^(−ΔCt). Bars show "
            "the housekeeping-normalised signal per target — *not* the raw "
            "ΔCt cycle counts. Use the camera icon on the figure toolbar "
            "to download any single plot as a PNG. To download every bar "
            "plot for one housekeeping gene as a single PNG, see the "
            "**Downloads** tab."
        ).classes("text-sm text-slate-600")

        for ref, res_df in results.items():
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
                _render_figures_for_results(
                    res_df,
                    value_col="Expr_vs_HK",
                    y_label="Relative expression (2^−ΔCt)",
                    ref=ref,
                    color_map=color_map,
                )


def _render_ddct_results(state: dict, refs: dict) -> None:
    panel = refs["ddct_panel"]
    panel.clear()
    results = state["ddct_results"] or {}
    ref_group = state["reference_group"]
    color_map = _color_map_for_results(state["standardized"], results)
    with panel:
        ui.markdown(
            f"**ΔΔCt** relative expression: 2^(−ΔΔCt), normalised so that "
            f"*{ref_group}* anchors at 1 within each batch. Use the camera "
            "icon on the figure toolbar to download any single plot as a "
            "PNG. To download every bar plot for one housekeeping gene as a "
            "single PNG, see the **Downloads** tab."
        ).classes("text-sm text-slate-600")

        for ref, res_df in results.items():
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
                _render_figures_for_results(
                    res_df,
                    value_col="Relative_Expr",
                    y_label="Relative expression (2^−ΔΔCt)",
                    ref=ref,
                    color_map=color_map,
                )


def _color_map_for_results(
    standardized: pd.DataFrame | None,
    results: dict[str, pd.DataFrame],
) -> dict[str, str]:
    """Build a single label→colour map shared across every plot.

    If the dataset has more than one biological group, colour by group;
    otherwise colour by sample. Items are taken from the standardised
    DataFrame so the order is deterministic and stable across plots.
    """
    groups: list[str] = []
    samples: list[str] = []
    if standardized is not None:
        groups = list(pd.unique(standardized["Group"].astype(str)))
        samples = list(pd.unique(standardized["Sample"].astype(str)))
    else:
        seen_g, seen_s = set(), set()
        for df in results.values():
            for g in df["Group"].astype(str):
                if g not in seen_g:
                    seen_g.add(g)
                    groups.append(g)
            for s in df["Sample"].astype(str):
                if s not in seen_s:
                    seen_s.add(s)
                    samples.append(s)
    items = groups if len(groups) > 1 else samples
    return _build_color_map(items)


def _render_figures_for_results(
    res_df: pd.DataFrame,
    value_col: str,
    y_label: str,
    ref: str,
    color_map: dict[str, str],
) -> None:
    targets = list(pd.unique(res_df["Target"]))
    groups = list(pd.unique(res_df["Group"]))
    use_group = len(groups) > 1
    cols = 1 if len(targets) == 1 else 2
    with ui.grid(columns=cols).classes("w-full gap-3"):
        for target in targets:
            sub = res_df[res_df["Target"] == target]
            fig = _plot_figure(
                sub, target, ref, use_group, value_col, y_label, color_map
            )
            ui.plotly(fig).classes("w-full h-80")


def _build_combined_figure_for_hk(
    ref: str,
    res_df: pd.DataFrame,
    value_col: str,
    y_label: str,
    color_map: dict[str, str],
) -> tuple[go.Figure, int]:
    """Combined subplot figure of every target for one housekeeping gene.

    Returns ``(figure, height_px)`` so the caller can size the plotly div
    sensibly. The user downloads it as a single PNG via the figure's
    built-in camera icon (rendered client-side by plotly.js — no Chrome
    needed server-side).
    """
    targets = list(pd.unique(res_df["Target"]))
    groups = list(pd.unique(res_df["Group"]))
    use_group = len(groups) > 1

    n = len(targets)
    cols = 1 if n == 1 else (2 if n <= 4 else 3)
    rows = max(1, (n + cols - 1) // cols)
    titles = [f"{t} / {ref}" for t in targets]
    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=titles,
        horizontal_spacing=0.08,
        vertical_spacing=0.18,
    )
    for i, target in enumerate(targets):
        r, c = i // cols + 1, i % cols + 1
        sub = res_df[res_df["Target"] == target]
        _add_bar_traces(fig, sub, use_group, value_col, color_map, row=r, col=c)
        fig.update_yaxes(title_text=y_label, row=r, col=c)

    height = max(320, 320 * rows)
    fig.update_layout(
        template="plotly_white",
        showlegend=False,
        plot_bgcolor="white",
        height=height,
        margin=dict(l=50, r=20, t=60, b=40),
        modebar=dict(remove=["lasso2d", "select2d"]),
        title=dict(
            text=f"{ref} · combined (download via camera icon)",
            x=0.5,
            xanchor="center",
        ),
    )
    return fig, height


def _render_downloads(state: dict, refs: dict) -> None:
    panel = refs["downloads_panel"]
    panel.clear()

    def _build_std() -> pd.DataFrame:
        std = state["standardized"].copy()
        std["Excluded"] = std["Well"].isin(state["excluded_wells"])
        # Honour the multi-batch checkbox: omit the Batch column from
        # exported sheets when the user said samples are *not* batched.
        if not state.get("has_batches") and "Batch" in std.columns:
            std = std.drop(columns=["Batch"])
        return std

    def _build_ddct() -> dict[str, pd.DataFrame]:
        # ddCt always carries a Batch column (defaulting to "batch_1") since
        # compute_delta_delta_ct is batch-aware. Strip it from the export
        # when the user said samples are *not* from different batches.
        results = state["ddct_results"] or {}
        if state.get("has_batches"):
            return results
        return {k: v.drop(columns=["Batch"], errors="ignore") for k, v in results.items()}

    with panel:
        ui.label("Excel workbook").classes("text-base font-semibold")
        ui.label(
            "Includes raw data (sorted), ΔCt and ΔΔCt sheets per "
            "housekeeping gene, plus the matching `formatted_*` sheets — "
            "wide grouped tables that show sample names alongside every value."
        ).classes("text-sm text-slate-600")

        def _download_xlsx() -> None:
            try:
                data = results_to_xlsx_bytes(
                    _build_std(),
                    state["dct_results"] or {},
                    _build_ddct(),
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"Excel build failed: {ex}", type="negative")
                return
            ui.download(data, f"{_filename_stem(state)}.xlsx")

        ui.button("Download xlsx", icon="download", on_click=_download_xlsx).props(
            "color=primary unelevated"
        )

        ui.separator()
        ui.label("CSV bundle (zip)").classes("text-base font-semibold")
        ui.label(
            "Same logical sheets as the Excel workbook, packed as one CSV "
            "per file inside a zip archive (raw_data.csv, dCt_*.csv, "
            "ddCt_*.csv, formatted_*.csv)."
        ).classes("text-sm text-slate-600")

        def _download_csv_zip() -> None:
            try:
                data = results_to_csv_zip_bytes(
                    _build_std(),
                    state["dct_results"] or {},
                    _build_ddct(),
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"CSV bundle build failed: {ex}", type="negative")
                return
            ui.download(data, f"{_filename_stem(state)}.zip")

        ui.button(
            "Download csv (zip)", icon="folder_zip", on_click=_download_csv_zip
        ).props("color=primary outline")

        ui.separator()
        ui.label("Figures").classes("text-base font-semibold")
        ui.label(
            "Each individual figure on the ΔCt / ΔΔCt result tabs has a "
            "built-in PNG download in its toolbar (camera icon). The "
            "combined views below pack every target for one housekeeping "
            "gene into a single PNG — click the camera icon on the "
            "respective figure to download."
        ).classes("text-sm text-slate-600")

        dct_results = state.get("dct_results") or {}
        ddct_results = state.get("ddct_results") or {}
        std = state.get("standardized")

        if dct_results:
            color_map = _color_map_for_results(std, dct_results)
            ui.label("Combined ΔCt plots").classes(
                "text-sm font-semibold mt-2"
            )
            for ref, res_df in dct_results.items():
                with ui.expansion(
                    f"ΔCt · {ref}", icon="show_chart"
                ).classes(
                    "w-full bg-slate-50 rounded border border-slate-200"
                ):
                    fig, height = _build_combined_figure_for_hk(
                        ref,
                        res_df,
                        "Expr_vs_HK",
                        "Relative expression (2^−ΔCt)",
                        color_map,
                    )
                    ui.label(
                        "Click the camera icon on the figure toolbar to "
                        "save every relative-expression bar plot for this "
                        "housekeeping gene as a single PNG."
                    ).classes("text-xs text-slate-600")
                    ui.plotly(fig).classes(
                        f"w-full h-[{height}px]"
                    )

        if ddct_results:
            color_map = _color_map_for_results(std, ddct_results)
            ui.label("Combined ΔΔCt plots").classes(
                "text-sm font-semibold mt-2"
            )
            for ref, res_df in ddct_results.items():
                with ui.expansion(
                    f"ΔΔCt · {ref}", icon="bar_chart"
                ).classes(
                    "w-full bg-slate-50 rounded border border-slate-200"
                ):
                    fig, height = _build_combined_figure_for_hk(
                        ref,
                        res_df,
                        "Relative_Expr",
                        "Relative expression (2^−ΔΔCt)",
                        color_map,
                    )
                    ui.label(
                        "Click the camera icon on the figure toolbar to "
                        "save every ΔΔCt bar plot for this housekeeping "
                        "gene as a single PNG."
                    ).classes("text-xs text-slate-600")
                    ui.plotly(fig).classes(
                        f"w-full h-[{height}px]"
                    )


def _add_bar_traces(
    fig: go.Figure,
    df: pd.DataFrame,
    use_group: bool,
    value_col: str,
    color_map: dict[str, str],
    *,
    row: int | None = None,
    col: int | None = None,
) -> None:
    """Add the (mean ± std bar) and optional sample-dot traces to ``fig``.

    Bars are coloured by their x-axis label using ``color_map``. The same
    map is used across every plot so a given group/sample always gets the
    same colour. Falls back to :data:`PRIMARY` for any label that was not
    pre-registered.
    """
    x_col = "Group" if use_group else "Sample"
    summary = (
        df.groupby(x_col, sort=False)[value_col]
        .agg(["mean", "std", "count"])
        .reset_index()
        .fillna(0)
    )
    bar_colors = [color_map.get(str(x), PRIMARY) for x in summary[x_col]]
    add_kw = {"row": row, "col": col} if row is not None else {}
    fig.add_bar(
        x=summary[x_col],
        y=summary["mean"],
        error_y=dict(
            type="data", array=summary["std"], visible=True, color="#334155"
        ),
        marker=dict(color=bar_colors, line=dict(color="#1e3a8a", width=0)),
        name="mean",
        hovertemplate="%{x}<br>mean = %{y:.3f}<extra></extra>",
        **add_kw,
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
            **add_kw,
        )


def _plot_figure(
    df: pd.DataFrame,
    target: str,
    ref: str,
    use_group: bool,
    value_col: str,
    y_label: str,
    color_map: dict[str, str],
) -> go.Figure:
    fig = go.Figure()
    _add_bar_traces(fig, df, use_group, value_col, color_map)
    fig.update_layout(
        title=dict(text=f"{target} / {ref}", x=0.5, xanchor="center"),
        yaxis_title=y_label,
        xaxis_title="",
        template="plotly_white",
        showlegend=False,
        margin=dict(l=40, r=10, t=50, b=40),
        plot_bgcolor="white",
        modebar=dict(remove=["lasso2d", "select2d"]),
    )
    return fig


def start(host: str | None = None, port: int | None = None) -> None:
    """Launch the NiceGUI server (blocking)."""
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
