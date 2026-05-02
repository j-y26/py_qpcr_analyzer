"""Desktop entry point — runs the NiceGUI app in a native pywebview window.

This is what the PyInstaller-bundled Windows .exe and macOS .app launch.
Unlike :mod:`qpcr_analyzer.__main__`, which serves on localhost for browser
access, this module opens a self-contained desktop window so non-technical
users never have to think about ports, browsers, or terminals.
"""

from __future__ import annotations

import multiprocessing
import os
import sys

from nicegui import ui

from qpcr_analyzer.app.main import index  # noqa: F401  -- registers the @ui.page("/")


def main() -> None:
    """Launch the app in a native desktop window (blocking)."""
    # NiceGUI picks a random free port internally when port=0 in native mode.
    ui.run(
        title="qPCR Analyzer",
        native=True,
        window_size=(1400, 900),
        reload=False,
        show=False,
        favicon="🧬",
        port=0,
    )


if __name__ in {"__main__", "__mp_main__"}:
    # PyInstaller-frozen apps must call freeze_support() before any
    # multiprocessing use (NiceGUI/uvicorn spawn workers on Windows).
    if getattr(sys, "frozen", False):
        multiprocessing.freeze_support()
        # Keep stdout/stderr writable when launched from a Windows .exe with
        # no console attached — uvicorn logs would otherwise crash on write.
        if sys.stdout is None:
            sys.stdout = open(os.devnull, "w")  # noqa: SIM115
        if sys.stderr is None:
            sys.stderr = open(os.devnull, "w")  # noqa: SIM115
    main()
