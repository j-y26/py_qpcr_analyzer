"""File-save helpers that work in both browser and native (pywebview) modes.

NiceGUI's ``ui.download`` triggers a browser download via a synthetic anchor.
That works in a real browser but is silently dropped by the embedded webview
runtimes (WKWebView on macOS, WebView2 on Windows) we use in ``desktop.py``,
because they don't have a "Downloads folder" concept the same way browsers do.

This module exposes :func:`save_bytes`, which prefers pywebview's native
"Save As…" dialog when running inside the bundled desktop app, and falls back
to ``ui.download`` for the standard browser-served deployment.

It also exposes :func:`install_native_plotly_bridge`, which monkey-patches
``Plotly.downloadImage`` on the client so the camera-icon button on every
Plotly figure routes back through the native save dialog instead of the
no-op blob download.
"""

from __future__ import annotations

import base64
from pathlib import Path

from nicegui import app, ui


def _native_window():
    """Return the pywebview window if running in native mode, else ``None``.

    NiceGUI exposes ``app.native.main_window`` only when ``ui.run(native=True)``
    is in effect; in browser-served mode the attribute access either raises or
    yields ``None``.
    """
    try:
        return app.native.main_window
    except (AttributeError, RuntimeError):
        return None


async def save_bytes(
    data: bytes,
    suggested_name: str,
    file_types: tuple[str, ...] = (),
) -> None:
    """Save *data* to disk, picking the right mechanism for the runtime.

    In native mode this opens pywebview's "Save As…" dialog and writes the
    bytes to the chosen path. In browser mode it falls back to
    ``ui.download``, which streams the bytes to the user's browser.

    *file_types* uses the pywebview format, e.g.
    ``("Excel workbook (*.xlsx)", "All files (*.*)")``.
    """
    window = _native_window()
    if window is None:
        ui.download(data, suggested_name)
        return

    # Imported lazily so the browser-only install (no pywebview) still works.
    import webview  # type: ignore[import-not-found]

    try:
        result = await window.create_file_dialog(
            dialog_type=webview.SAVE_DIALOG,
            save_filename=suggested_name,
            file_types=file_types,
        )
    except Exception as ex:  # noqa: BLE001 -- surface any dialog error to UI
        ui.notify(f"Save dialog failed: {ex}", type="negative")
        return

    if not result:
        return  # user cancelled

    path = result if isinstance(result, str) else result[0]
    try:
        Path(path).write_bytes(data)
    except OSError as ex:
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
    if _native_window() is None:
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
