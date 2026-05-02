"""qPCR Analyzer — lightweight browser-based qPCR data analysis.

Public surface
--------------
- :mod:`qpcr_analyzer.core`     pure-Python analysis library (no UI deps).
- :mod:`qpcr_analyzer.app.main` NiceGUI front-end; exposes ``start(host, port)``.
- :mod:`qpcr_analyzer.__main__` entry point bound to the ``qpcr-analyzer``
  console script declared in ``pyproject.toml``.

The package is split so that ``core`` can be imported and used as a library
without pulling in any UI dependencies.
"""

__version__ = "1.1.0"

