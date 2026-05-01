"""NiceGUI front-end for qPCR Analyzer.

The page is a two-pane layout:

* **Left pane** — a vertical stepper that drives the workflow:
  upload → column mapping → groups & batches → outliers & per-HK exclusion
  → run ΔCt → reference group & run ΔΔCt.
* **Right pane** — persistent output tabs that activate as data flows in:
  *Summary*, *Data preview*, *Excluded blocks*, *ΔCt results*, *ΔΔCt
  results*, *Downloads*.

All scientific logic lives in :mod:`qpcr_analyzer.core`; this module only
collects user input, calls the core, and renders Plotly figures + tables.
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
    build_blocks,
    compute_delta_ct,
    compute_delta_delta_ct,
    detect_columns,
    mark_outliers,
    read_table,
    results_to_csv_zip_bytes,
    results_to_xlsx_bytes,
    samples_missing_hk,
    sort_wells,
    summarize_dataset,
    target_order,
    validate_sample_batches,
    validate_sample_groups,
)
from qpcr_analyzer.core.columns import REQUIRED

PRIMARY = "#2563eb"
PALETTE = ["#2563eb", "#0ea5e9", "#14b8a6", "#f59e0b", "#ef4444", "#a855f7"]


def _new_state() -> dict:
    """Build a fresh per-session state dict."""
    return {
        "filename": None,
        "raw_df": None,
        "mapping": None,
        "standardized": None,
        "summary": None,
        "sample_batches": {},        # {sample: batch_label}
        "reference_group": None,
        "tolerance": 1.0,
        "flagged": None,
        "excluded_wells": set(),
        "_outliers_initialized": False,
        "ref_genes": [],
        # {hk_gene: set(sample_name, ...)}, "*" means all HKs
        "sample_excludes_per_hk": {},
        "dct_results": None,
        "ddct_results": None,
        "dct_done": False,
    }


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
            ui.label("· ΔCt & batch-aware ΔΔCt quantification").classes(
                "text-sm text-slate-500"
            )
        ui.space()
        dark = ui.dark_mode()
        ui.button(icon="dark_mode", on_click=dark.toggle).props("flat round dense")

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
                    "4. **Mark outliers and pick housekeeping gene(s).** "
                    "Samples missing a valid HK Cq are surfaced for confirmation "
                    "before being skipped for that HK only.\n"
                    "5. **Run ΔCt.** Always required.\n"
                    "6. **Run ΔΔCt** (optional). Pick a reference group and "
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
            # Seed sample_batches from the file if a Batch column is present.
            if "Batch" in std.columns:
                state["sample_batches"] = (
                    std.groupby("Sample")["Batch"].first().astype(str).to_dict()
                )
            else:
                state["sample_batches"] = {}
            state["reference_group"] = None
            state["dct_done"] = False
            state["dct_results"] = None
            state["ddct_results"] = None
            state["summary"] = summarize_dataset(std, filename=state["filename"])
            _render_full_summary(state, refs)
            _render_data_preview(state, refs)
            _render_groups(state, refs)
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            ui.button("Continue", on_click=go_to_groups).props(
                "color=primary unelevated"
            )


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
            # Collect inline edits from the table
            edited = refs.get("groups_rows") or []
            if edited:
                rows_df = pd.DataFrame(edited)
                std = state["standardized"].copy()
                # Apply group / batch overrides keyed on Sample
                group_map = dict(zip(rows_df["Sample"], rows_df["Group"]))
                batch_map = dict(zip(rows_df["Sample"], rows_df["Batch"]))
                std["Group"] = std["Sample"].map(group_map).fillna(std["Group"]).astype(str)
                std["Batch"] = std["Sample"].map(batch_map).fillna("batch_1").astype(str)
                state["standardized"] = std
                state["sample_batches"] = {
                    str(s): str(b) for s, b in batch_map.items()
                }
                state["summary"] = summarize_dataset(std, filename=state["filename"])
                _render_full_summary(state, refs)
                _render_data_preview(state, refs)

            sg_errs = validate_sample_groups(state["standardized"])
            sb_errs = validate_sample_batches(state["standardized"])
            for msg in (*sg_errs, *sb_errs):
                ui.notify(msg, type="negative")
            if sg_errs or sb_errs:
                return

            _render_outliers(state, refs)
            stepper.next()

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            ui.button("Continue", on_click=go_to_outliers).props(
                "color=primary unelevated"
            )


def _build_step_outliers(state: dict, refs: dict, stepper) -> None:
    with ui.step(
        "outliers",
        title="4. Outliers & housekeeping gene(s)",
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

        with ui.stepper_navigation():
            ui.button("Back", on_click=stepper.previous).props("flat")
            ui.button("Continue", on_click=stepper.next).props(
                "color=primary unelevated"
            )


def _build_step_dct(state: dict, refs: dict, stepper) -> None:
    with ui.step("dct", title="5. Run ΔCt", icon="play_arrow"):
        ui.markdown(
            "ΔCt normalises each sample to its housekeeping gene Cq. "
            "No reference group is required at this step — the result is "
            "ready as soon as you click below.\n\n"
            "Samples flagged in the panel above as *missing HK Cq* will be "
            "skipped for that HK gene only. Use the per-HK exclusion controls "
            "to drop additional samples."
        ).classes("text-slate-700")
        refs["dct_status"] = ui.label("").classes("text-sm text-slate-600")

        def run_dct() -> None:
            if state["standardized"] is None:
                ui.notify("Load and configure data first.", type="negative")
                return
            if not state.get("ref_genes"):
                ui.notify("Pick at least one housekeeping gene.", type="negative")
                return
            std = state["standardized"].copy()
            std["Excluded"] = std["Well"].isin(state["excluded_wells"])
            try:
                dct = compute_delta_ct(
                    std,
                    list(state["ref_genes"]),
                    sample_excludes_per_hk=state["sample_excludes_per_hk"],
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"ΔCt failed: {ex}", type="negative")
                return
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
            ui.button("Run ΔCt", on_click=run_dct).props(
                "color=primary unelevated icon=play_arrow"
            )


def _build_step_ddct(state: dict, refs: dict, stepper) -> None:
    with ui.step("ddct", title="6. Run ΔΔCt (optional)", icon="bar_chart"):
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

        def run_ddct() -> None:
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
            try:
                ddct = compute_delta_delta_ct(
                    std,
                    list(state["ref_genes"]),
                    reference_group=ref_group,
                    sample_batches=state["sample_batches"] or None,
                    sample_excludes_per_hk=state["sample_excludes_per_hk"],
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"ΔΔCt failed: {ex}", type="negative")
                return
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
        ui.label("ΔΔCt configuration").classes("text-base font-semibold")
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
    with panel:
        ui.label("Dataset summary").classes("text-base font-semibold")
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
        _summary_section("Samples", s["n_samples"], s["samples"])
        _summary_section("Targets", s["n_targets"], s["targets"])
        _summary_section("Groups", s["n_groups"], s["groups"])
        if s["has_batch_column"]:
            _summary_section("Batches (from file)", s["n_batches"], s["batches"])
        else:
            ui.label(
                "No Batch column in file — assign batches in step 3 if needed."
            ).classes("text-xs text-slate-500")


def _stat_chip(label: str, value, warn: bool = False) -> None:
    color = "amber" if warn else "blue"
    with ui.element("div").classes(
        f"px-3 py-1 rounded-full bg-{color}-50 border border-{color}-200 "
        "text-sm flex items-center gap-1"
    ):
        ui.label(str(value)).classes(f"font-semibold text-{color}-700")
        ui.label(label).classes("text-slate-600")


def _summary_section(title: str, count: int, values: list[str]) -> None:
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
    cols = [c for c in ["Well", "Target", "Group", "Sample", "Batch", "Cq"] if c in sorted_df.columns]
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

    has_batch = "Batch" in df.columns
    if has_batch:
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
    # Initial dropdown options — combobox lets users either pick from these
    # or type a new value, which then becomes available to subsequent rows.
    initial_groups = sorted({str(r["Group"]) for r in rows})
    initial_batches = sorted({str(r["Batch"]) for r in rows})

    with container:
        ui.label(
            f"{len(samples_meta)} unique sample(s). "
            f"{'Batch column detected — pre-filled from file.' if has_batch else 'No Batch column in file.'}"
        ).classes("text-sm text-slate-600")
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

        ui.label(
            "Tip: leave Batch as 'batch_1' for a single experimental run."
        ).classes("text-xs text-slate-400")


def _render_outliers(state: dict, refs: dict) -> None:
    container = refs["outlier_container"]
    container.clear()
    df = state["standardized"]
    targets = target_order(df)

    with container:
        with ui.row().classes("w-full items-end gap-4 flex-wrap"):
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
                _rebuild_outlier_view(state, refs)

            ui.button("Apply", on_click=_apply).props("color=primary outline")

        refs["outlier_summary_slot"] = ui.label("").classes("text-sm text-slate-600")
        refs["excluded_select_slot"] = ui.column().classes("w-full")
        refs["per_hk_slot"] = ui.column().classes("w-full")

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
            _refresh_per_hk_panel(state, refs)

        excluded_select.on_value_change(_on_change)

    refs["outlier_summary_slot"].set_text(_excluded_summary(state))
    _render_excluded_blocks(state, refs)
    _refresh_per_hk_panel(state, refs)


def _refresh_per_hk_panel(state: dict, refs: dict) -> None:
    """Render the per-HK sample exclusion controls below the main outlier UI.

    For each housekeeping gene the user has selected we list samples that
    have no usable HK Cq (auto-flagged for exclusion) and provide a
    multi-select for additional manual sample exclusions, plus a separate
    "exclude entirely (all HKs)" multi-select.
    """
    slot = refs["per_hk_slot"]
    slot.clear()
    if not state.get("ref_genes"):
        with slot:
            ui.label(
                "Pick a housekeeping gene above to configure per-HK sample "
                "exclusion."
            ).classes("text-xs text-slate-500")
        return

    std = state["standardized"].copy()
    std["Excluded"] = std["Well"].isin(state["excluded_wells"])
    missing = samples_missing_hk(std, list(state["ref_genes"]))
    all_samples = sorted(std["Sample"].astype(str).unique())

    excludes = state.setdefault("sample_excludes_per_hk", {})

    with slot:
        ui.label("Per-HK sample exclusion").classes(
            "text-base font-semibold mt-2"
        )
        with ui.expansion("Why exclude samples per HK?", icon="info").classes(
            "w-full"
        ):
            ui.markdown(
                "ΔCt is computed *per housekeeping gene*. If a sample has no "
                "valid HK Cq for one gene, it can still be analysed against "
                "another HK. Excluding a sample for a single HK keeps it in "
                "the rest of the analysis instead of breaking the whole run."
            ).classes("text-sm text-slate-600")

        # Auto-fill missing-HK samples and surface them for confirmation
        for hk in state["ref_genes"]:
            auto_excl = set(missing.get(hk, []))
            current = set(excludes.get(hk, set())) | auto_excl
            excludes[hk] = current

            with ui.card().classes("w-full border border-slate-200 shadow-none"):
                ui.label(f"Housekeeping: {hk}").classes("font-semibold")
                if auto_excl:
                    ui.label(
                        "Auto-excluded (no valid HK Cq): "
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

                def _on_change(_e, hk=hk, sel=sel) -> None:
                    excludes[hk] = set(sel.value or [])

                sel.on_value_change(_on_change)

        # Global "exclude entirely" control
        with ui.card().classes("w-full border border-slate-200 shadow-none"):
            ui.label("Exclude entirely (all housekeeping genes)").classes(
                "font-semibold"
            )
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
                excludes["*"] = set(sel.value or [])

            global_sel.on_value_change(_on_global)


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
    with panel:
        ui.markdown(
            "**ΔCt** = mean Cq(target) − mean Cq(housekeeping gene), per "
            "sample. Use the camera icon on the figure toolbar to download "
            "any plot as a PNG."
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
                    res_df, value_col="dCt", y_label="ΔCt", ref=ref
                )


def _render_ddct_results(state: dict, refs: dict) -> None:
    panel = refs["ddct_panel"]
    panel.clear()
    results = state["ddct_results"] or {}
    ref_group = state["reference_group"]
    with panel:
        ui.markdown(
            f"**ΔΔCt** relative expression: 2^(−ΔΔCt), normalised so that "
            f"*{ref_group}* anchors at 1 within each batch. Use the camera "
            "icon on the figure toolbar to download any plot as a PNG."
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
                )


def _render_figures_for_results(
    res_df: pd.DataFrame, value_col: str, y_label: str, ref: str
) -> None:
    targets = list(pd.unique(res_df["Target"]))
    groups = list(pd.unique(res_df["Group"]))
    use_group = len(groups) > 1
    cols = 1 if len(targets) == 1 else 2
    with ui.grid(columns=cols).classes("w-full gap-3"):
        for target in targets:
            sub = res_df[res_df["Target"] == target]
            fig = _plot_figure(sub, target, ref, use_group, value_col, y_label)
            ui.plotly(fig).classes("w-full h-80")


def _render_downloads(state: dict, refs: dict) -> None:
    panel = refs["downloads_panel"]
    panel.clear()

    def _build_std() -> pd.DataFrame:
        std = state["standardized"].copy()
        std["Excluded"] = std["Well"].isin(state["excluded_wells"])
        return std

    with panel:
        ui.label("Excel workbook").classes("text-base font-semibold")
        ui.label(
            "Includes raw data (sorted), ΔCt and ΔΔCt sheets per HK, plus "
            "the matching `formatted_*` sheets — wide grouped tables that "
            "show sample names alongside every value."
        ).classes("text-sm text-slate-600")

        def _download_xlsx() -> None:
            try:
                data = results_to_xlsx_bytes(
                    _build_std(),
                    state["dct_results"] or {},
                    state["ddct_results"] or {},
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"Excel build failed: {ex}", type="negative")
                return
            stem = Path(state["filename"] or "data").stem
            ui.download(data, f"qpcr_results_{stem}.xlsx")

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
                    state["ddct_results"] or {},
                )
            except Exception as ex:  # noqa: BLE001
                ui.notify(f"CSV bundle build failed: {ex}", type="negative")
                return
            stem = Path(state["filename"] or "data").stem
            ui.download(data, f"qpcr_results_{stem}.zip")

        ui.button(
            "Download csv (zip)", icon="folder_zip", on_click=_download_csv_zip
        ).props("color=primary outline")

        ui.separator()
        ui.label("Figures").classes("text-base font-semibold")
        ui.label(
            "Each figure on the ΔCt and ΔΔCt result tabs has a built-in PNG "
            "download in its toolbar (camera icon)."
        ).classes("text-sm text-slate-600")


def _plot_figure(
    df: pd.DataFrame,
    target: str,
    ref: str,
    use_group: bool,
    value_col: str,
    y_label: str,
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
        hovertemplate="%{x}<br>mean = %{y:.3f}<extra></extra>",
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
        xaxis_title="",
        template="plotly_white",
        showlegend=False,
        margin=dict(l=40, r=10, t=50, b=40),
        plot_bgcolor="white",
    )
    fig.update_layout(
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
