"""Dev-only launcher with hot-reload. Not shipped in the package."""
from nicegui import ui
from qpcr_analyzer.app.main import index  # noqa: F401  -- registers the @ui.page

ui.run(host="127.0.0.1", port=8090, reload=True, show=False, title="qPCR Analyzer (dev)")
