"""File-save helpers that work in both browser and native (pywebview) modes.

NiceGUI's ``ui.download`` triggers a browser download via a synthetic anchor.
That works in a real browser but is silently dropped by the embedded webview
runtimes (WKWebView on macOS, WebView2 on Windows) we use in ``desktop.py``,
because they don't have a "Downloads folder" concept the same way browsers do.

This module exposes :func:`save_bytes`, which prefers a native ``tkinter``
"Save As…" dialog when running inside the bundled desktop app, and falls back
to ``ui.download`` for the standard browser-served deployment.

Why ``tkinter`` and not ``pywebview.create_file_dialog``?
The pywebview dialog is a synchronous call that must run on the platform GUI
thread; NiceGUI's async wrapper around it has been flaky in practice (silent
no-op on some webview backends). ``tkinter.filedialog`` is stdlib, runs on a
private GUI thread we control, and behaves identically on Windows and macOS.

It also exposes :func:`install_native_plotly_bridge`, which monkey-patches
``Plotly.downloadImage`` on the client so the camera-icon button on every
Plotly figure routes back through the native save dialog instead of the
no-op blob download.
"""

from __future__ import annotations

import asyncio
import base64
import logging
import os
import sys
from pathlib import Path

from nicegui import app, ui

logger = logging.getLogger(__name__)


def _running_native() -> bool:
    """True when running inside the bundled desktop app or with ``native=True``.

    ``sys.frozen`` is set by PyInstaller in the shipped .exe / .app bundles,
    so we can decide without depending on NiceGUI's ``app.native.main_window``,
    which has timing-dependent semantics.
    """
    if getattr(sys, "frozen", False):
        return True
    try:
        return app.native.main_window is not None
    except (AttributeError, RuntimeError):
        return False


def _ask_save_path(
    suggested_name: str,
    tk_file_types: list[tuple[str, str]],
) -> str | None:
    """Show a tkinter Save-As dialog. Blocking — call via ``asyncio.to_thread``."""
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    # Force the dialog above the webview window — Windows/macOS sometimes
    # park new top-levels behind the focused webview otherwise.
    try:
        root.attributes("-topmost", True)
    except tk.TclError:
        pass
    try:
        ext = os.path.splitext(suggested_name)[1]
        path = filedialog.asksaveasfilename(
            parent=root,
            initialfile=suggested_name,
            filetypes=tk_file_types or [("All files", "*.*")],
            defaultextension=ext,
        )
    finally:
        root.update_idletasks()
        root.destroy()
    return path or None


def _normalize_file_types(specs: tuple[str, ...]) -> list[tuple[str, str]]:
    """Translate ``"Label (*.ext)"`` strings into tk's ``(label, "*.ext")``."""
    out: list[tuple[str, str]] = []
    for spec in specs:
        if "(" in spec and ")" in spec:
            label, _, rest = spec.partition("(")
            ext = rest.rstrip(")").strip()
            out.append((label.strip(), ext or "*.*"))
        else:
            out.append((spec, "*.*"))
    return out


async def save_bytes(
    data: bytes,
    suggested_name: str,
    file_types: tuple[str, ...] = (),
) -> None:
    """Save *data* to disk, picking the right mechanism for the runtime.

    In native mode this opens a ``tkinter`` "Save As…" dialog and writes the
    bytes to the chosen path. In browser mode it falls back to
    ``ui.download``, which streams the bytes to the user's browser.

    *file_types* uses ``"Label (*.ext)"`` strings, e.g.
    ``("Excel workbook (*.xlsx)", "All files (*.*)")``.
    """
    if not _running_native():
        ui.download(data, suggested_name)
        return

    tk_types = _normalize_file_types(file_types)
    try:
        path = await asyncio.to_thread(_ask_save_path, suggested_name, tk_types)
    except Exception as ex:  # noqa: BLE001 -- surface any dialog error to UI
        logger.exception("save-as dialog failed")
        ui.notify(f"Save dialog failed: {ex}", type="negative")
        return

    if not path:
        return  # user cancelled

    try:
        Path(path).write_bytes(data)
    except OSError as ex:
        logger.exception("file write failed")
        ui.notify(f"Could not write {path}: {ex}", type="negative")
        return

    ui.notify(f"Saved to {path}", type="positive")


# JS shim that overrides Plotly.downloadImage so the modebar camera icon
# round-trips through Python instead of attempting a webview download.
# Plotly may load after our script (NiceGUI lazy-imports it the first time
# a figure renders), so we poll for the global until it's defined.
_PLOTLY_PATCH_JS = """
(function() {
  function patch() {
    if (!window.Plotly) { setTimeout(patch, 100); return; }
    if (window.Plotly._qpcrPatched) return;
    window.Plotly._qpcrPatched = true;
    const origDownload = window.Plotly.downloadImage;
    window.Plotly.downloadImage = async function(gd, opts) {
      opts = opts || {};
      const fmt = opts.format || 'png';
      try {
        const dataUrl = await window.Plotly.toImage(gd, {
          format: fmt,
          width:  opts.width  || (gd._fullLayout && gd._fullLayout.width)  || 1024,
          height: opts.height || (gd._fullLayout && gd._fullLayout.height) || 768,
          scale:  opts.scale  || 2,
        });
        const name = (opts.filename || 'plot') + '.' + fmt;
        emitEvent('qpcr-save-image', {name: name, dataUrl: dataUrl});
      } catch (e) {
        console.error('qPCR Plotly download patch failed', e);
        return origDownload.call(this, gd, opts);
      }
    };
  }
  patch();
})();
"""


def install_native_plotly_bridge() -> None:
    """Wire up the client→server save bridge for Plotly's camera button.

    Called from the page handler. Idempotent — only acts in native mode.
    The JS shim emits ``qpcr-save-image`` with ``{name, dataUrl}``; the
    Python handler decodes the base64 payload and routes through
    :func:`save_bytes`.
    """
    if not _running_native():
        return

    ui.add_body_html(f"<script>{_PLOTLY_PATCH_JS}</script>")

    async def _handle(event) -> None:
        # NiceGUI delivers the JS payload as event.args (dict for object args).
        args = event.args if hasattr(event, "args") else event
        name = (args or {}).get("name") or "plot.png"
        data_url = (args or {}).get("dataUrl") or ""
        if "," not in data_url:
            ui.notify("Plotly returned an empty image", type="negative")
            return
        _, b64 = data_url.split(",", 1)
        try:
            data = base64.b64decode(b64)
        except (ValueError, TypeError) as ex:
            ui.notify(f"Could not decode image: {ex}", type="negative")
            return
        ext = name.rsplit(".", 1)[-1].upper() if "." in name else "PNG"
        await save_bytes(
            data,
            name,
            file_types=(f"{ext} image (*.{ext.lower()})", "All files (*.*)"),
        )

    ui.on("qpcr-save-image", _handle)
