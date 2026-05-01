"""Console-script entry point for ``qpcr-analyzer``.

Invoked in three equivalent ways::

    qpcr-analyzer              # console script declared in pyproject.toml
    python -m qpcr_analyzer    # module-style invocation
    python src/qpcr_analyzer/__main__.py

All three call :func:`qpcr_analyzer.app.main.start`, which reads the
``QPCR_HOST`` / ``QPCR_PORT`` environment variables and launches NiceGUI.

The dual ``__main__`` / ``__mp_main__`` guard is required because NiceGUI's
hot-reload and Windows multiprocessing both re-import the entry module under
the alternate name ``__mp_main__``.
"""

from qpcr_analyzer.app.main import start


def main() -> None:
    """Launch the NiceGUI server. See :func:`qpcr_analyzer.app.main.start`."""
    start()


if __name__ in {"__main__", "__mp_main__"}:
    main()

