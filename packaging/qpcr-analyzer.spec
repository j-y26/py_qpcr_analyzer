# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for qPCR Analyzer desktop builds.

Produces a one-folder build because NiceGUI ships static assets (HTML
templates, JS, fonts) that one-file mode unpacks to a temp dir on every
launch — slow startup and antivirus false-positives. One-folder is
faster, easier to sign, and what the NiceGUI docs recommend.

Invoke from the repo root:

    pyinstaller packaging/qpcr-analyzer.spec --noconfirm

The CI workflow then wraps the resulting ``dist/qpcr-analyzer/`` folder
in an Inno Setup installer (Windows) or a .dmg (macOS).
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

# Spec files don't expose __file__; PyInstaller injects SPECPATH instead.
spec_root = Path(SPECPATH).resolve()  # noqa: F821 -- SPECPATH is injected
project_root = spec_root.parent
entry = project_root / "src" / "qpcr_analyzer" / "desktop.py"

# NiceGUI ships HTML/JS/CSS as data files that PyInstaller can't auto-detect
# from import analysis. Same story for plotly's JSON schema and pywebview's
# platform-specific runtime glue.
datas = []
datas += collect_data_files("nicegui")
datas += collect_data_files("plotly")
datas += collect_data_files("webview")  # pywebview package name
datas += collect_data_files("openpyxl")

# Some NiceGUI/uvicorn dynamic imports aren't picked up by static analysis.
hiddenimports = []
hiddenimports += collect_submodules("nicegui")
hiddenimports += collect_submodules("uvicorn")
hiddenimports += collect_submodules("webview")
hiddenimports += [
    "qpcr_analyzer",
    "qpcr_analyzer.app.main",
    "qpcr_analyzer.app._save",
    "qpcr_analyzer.core",
    "qpcr_analyzer.core.columns",
    "qpcr_analyzer.core.export",
    "qpcr_analyzer.core.io",
    "qpcr_analyzer.core.outliers",
    "qpcr_analyzer.core.quant",
    "qpcr_analyzer.core.summary",
    # tkinter ships with CPython but PyInstaller's analysis sometimes misses
    # filedialog when it's only referenced from a function imported at use
    # site — list explicitly so the Save-As dialog works in the bundle.
    "tkinter",
    "tkinter.filedialog",
]

block_cipher = None


a = Analysis(  # noqa: F821
    [str(entry)],
    pathex=[str(project_root / "src")],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Trim heavy unused libs that pandas/numpy may transitively pull in.
    # NB: tkinter is intentionally NOT excluded — the desktop "Save As…"
    # dialog (qpcr_analyzer.app._save) uses tkinter.filedialog.
    excludes=["PyQt5", "PyQt6", "PySide2", "PySide6", "matplotlib", "scipy", "IPython"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="qpcr-analyzer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,           # UPX often trips antivirus heuristics; skip it.
    console=False,       # GUI app — no console window on Windows.
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,    # CI sets this per-runner via --target-arch flag.
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)

coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="qpcr-analyzer",
)

# macOS .app bundle. Ignored on Windows/Linux builds.
app = BUNDLE(  # noqa: F821
    coll,
    name="qPCR Analyzer.app",
    icon=None,
    bundle_identifier="ca.sickkids.qpcr-analyzer",
    info_plist={
        "CFBundleName": "qPCR Analyzer",
        "CFBundleDisplayName": "qPCR Analyzer",
        "CFBundleShortVersionString": "2.1.2",
        "CFBundleVersion": "2.1.2",
        "NSHighResolutionCapable": True,
        # Required so pywebview's WKWebView can talk to the bundled localhost
        # NiceGUI server without ATS blocking it.
        "NSAppTransportSecurity": {
            "NSAllowsLocalNetworking": True,
        },
        "LSMinimumSystemVersion": "10.15.0",
    },
)
