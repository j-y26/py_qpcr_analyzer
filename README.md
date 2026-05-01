# qPCR Analyzer

A lightweight, browser-based tool for relative quantification of RT-qPCR data.
Supports both **ΔCt** (housekeeping-gene normalisation) and **batch-aware ΔΔCt**
(relative expression against a reference biological group).

* Pure-Python core (`pandas` / `numpy`) with a [NiceGUI](https://nicegui.io/)
  front-end — runs locally as a small web app, no internet connection needed.
* Cross-platform: Windows, macOS, and Linux. Anywhere CPython 3.10+ runs.
* No database, no cloud, no telemetry. Your data stays on your machine.

---

## Table of contents

1. [Quick start](#quick-start)
2. [Installation](#installation)
   - [Windows](#windows)
   - [macOS](#macos)
   - [Linux](#linux)
3. [Running the application](#running-the-application)
   - [Test environment](#test-environment-developmenttrial)
   - [Production environment](#production-environment-lab-server)
4. [Usage walkthrough](#usage-walkthrough)
5. [Methods](#methods)
6. [Output workbook](#output-workbook)
7. [Programmatic API](#programmatic-api)
8. [Project layout](#project-layout)
9. [Testing](#testing)
10. [License](#license)

---

## Quick start

Once installed (see below), run:

```bash
qpcr-analyzer
```

then open <http://127.0.0.1:8090> in your browser.

---

## Installation

### Prerequisites

| Requirement | Why |
|-------------|-----|
| Python ≥ 3.10 | minimum supported runtime |
| `pip` ≥ 23 | for `pyproject.toml` builds |
| (optional) `git` | only needed if installing from source |

Check with:

```bash
python --version
pip --version
```

> **Tip — virtual environments.** We strongly recommend installing into a
> virtual environment so qPCR Analyzer's dependencies don't collide with
> other Python tools on your system. Commands below assume one is active.

### Windows

PowerShell, from any folder you like:

```powershell
# 1. Install Python 3.10+ from https://www.python.org/downloads/windows/
#    (tick "Add python.exe to PATH" in the installer)

# 2. Create and activate a virtual environment
python -m venv qpcr-venv
qpcr-venv\Scripts\Activate.ps1

# 3. Install qPCR Analyzer from GitHub
pip install "git+https://github.com/<your-org>/qpcr-analyzer.git"

# 4. Run it
qpcr-analyzer
```

If PowerShell blocks the activation script, run once:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

### macOS

```bash
# 1. Install Python (Homebrew is easiest)
brew install python@3.12

# 2. Create and activate a virtual environment
python3 -m venv qpcr-venv
source qpcr-venv/bin/activate

# 3. Install qPCR Analyzer from GitHub
pip install "git+https://github.com/<your-org>/qpcr-analyzer.git"

# 4. Run it
qpcr-analyzer
```

On Apple Silicon, all required wheels (`pandas`, `numpy`, `openpyxl`,
`plotly`, `rapidfuzz`, `nicegui`) are published as native arm64, so
installation does not need to compile from source.

### Linux

```bash
# 1. Make sure python3-venv is available (Debian/Ubuntu)
sudo apt install python3 python3-venv python3-pip      # Debian/Ubuntu
# sudo dnf install python3 python3-pip                  # Fedora/RHEL
# sudo pacman -S python python-pip                      # Arch

# 2. Create and activate a virtual environment
python3 -m venv qpcr-venv
source qpcr-venv/bin/activate

# 3. Install qPCR Analyzer from GitHub
pip install "git+https://github.com/<your-org>/qpcr-analyzer.git"

# 4. Run it
qpcr-analyzer
```

### Optional extras

| Extra  | Purpose | Install command |
|--------|---------|-----------------|
| `xls`  | Read legacy `.xls` files (Excel 97-2003) | `pip install "qpcr-analyzer[xls] @ git+https://github.com/<your-org>/qpcr-analyzer.git"` |
| `dev`  | Test runner + linter for contributors | `pip install -e ".[dev]"` (after cloning) |

### Installing a specific version

```bash
pip install "git+https://github.com/<your-org>/qpcr-analyzer.git@v1.0.0"
```

### Upgrading

```bash
pip install --upgrade --force-reinstall \
    "git+https://github.com/<your-org>/qpcr-analyzer.git"
```

### Uninstalling

```bash
pip uninstall qpcr-analyzer
```

---

## Running the application

### Test environment (development/trial)

Use this when you want to try the app, run the test suite, or modify the
source code.

```bash
git clone https://github.com/<your-org>/qpcr-analyzer.git
cd qpcr-analyzer

# create + activate a venv (see platform sections above)
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\Activate.ps1

# editable install with dev extras
pip install -e ".[dev]"

# run the unit tests
pytest

# launch the app (default: http://127.0.0.1:8090, loopback only)
qpcr-analyzer
```

Editable mode (`-e`) means your changes to `src/qpcr_analyzer/` take effect
the next time you restart the server, no reinstall needed.

### Production environment (lab server)

Use this when you want a stable instance that colleagues on your LAN can
access. We recommend running it as a managed service so it restarts
automatically.

```bash
# 1. Install into a clean venv (system-wide is fine but not required)
python3 -m venv /opt/qpcr-analyzer
/opt/qpcr-analyzer/bin/pip install \
    "qpcr-analyzer[xls] @ git+https://github.com/<your-org>/qpcr-analyzer.git"

# 2. Run the server bound to all interfaces on the standard port
QPCR_HOST=0.0.0.0 QPCR_PORT=8090 /opt/qpcr-analyzer/bin/qpcr-analyzer
```

Then point colleagues at `http://<server-hostname>:8090`.

#### Configuration (environment variables)

| Variable    | Default     | Description                                          |
|-------------|-------------|------------------------------------------------------|
| `QPCR_HOST` | `127.0.0.1` | Bind address. Use `0.0.0.0` to expose on the LAN.    |
| `QPCR_PORT` | `8090`      | TCP port. (Default avoids clashing with Thermo Fisher's *Design & Analysis* software, which listens on 8080.) |

#### Run as a systemd service (Linux)

Create `/etc/systemd/system/qpcr-analyzer.service`:

```ini
[Unit]
Description=qPCR Analyzer
After=network.target

[Service]
Type=simple
Environment=QPCR_HOST=0.0.0.0
Environment=QPCR_PORT=8090
ExecStart=/opt/qpcr-analyzer/bin/qpcr-analyzer
Restart=on-failure
User=qpcr

[Install]
WantedBy=multi-user.target
```

Then:

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now qpcr-analyzer
sudo systemctl status qpcr-analyzer
```

#### Run as a Windows service

Use [NSSM](https://nssm.cc/) — install it, then:

```powershell
nssm install qpcr-analyzer "C:\opt\qpcr-analyzer\Scripts\qpcr-analyzer.exe"
nssm set qpcr-analyzer AppEnvironmentExtra QPCR_HOST=0.0.0.0 QPCR_PORT=8090
nssm start qpcr-analyzer
```

#### Hardening notes

* The server has **no built-in authentication**. Run it on a trusted
  internal network or front it with an authenticating reverse proxy
  (nginx + `auth_basic`, Caddy + Authelia, …).
* Use `https` via your reverse proxy if exposing beyond the host machine.
* Every uploaded file stays in memory only; nothing is written to disk
  unless the user clicks **Download xlsx**.

---

## Usage walkthrough

The interface is a five-step stepper:

| Step | What you do |
|------|-------------|
| **1. Upload**             | Drop in an `.xlsx`, `.xls`, `.csv`, `.tsv`, or `.txt` file (one row per well). |
| **2. Column mapping**     | Confirm or adjust auto-detected role assignments (Well, Target, Sample, Cq, Group). |
| **3. Groups & batches**   | Review which samples belong to which biological group. Optionally assign samples to **batches** if they come from multiple independent runs. |
| **4. Reference & outliers** | Choose the **reference group** (ΔΔCt anchor), **housekeeping gene(s)**, and replicate tolerance. Exclude outlier wells. |
| **5. Results**            | Inspect ΔCt and ΔΔCt bar charts. Download the Excel workbook. |

### Column auto-detection

The detector scores every column name against a synonym list using exact,
substring, and fuzzy matching (`rapidfuzz`). Columns scoring below 0.85
are left unassigned and shown as "unmatched" in the UI.

| Role   | Example column names                              |
|--------|---------------------------------------------------|
| Well   | Well, Well Position, Location                     |
| Target | Target, Gene, Assay, Detector                     |
| Sample | Sample, Sample ID, Sample Name                    |
| Cq     | Cq, Ct, Cq Value, Ct Mean                         |
| Group  | Group, Condition, Treatment, Biological Set Name  |

**Applied Biosystems note:** files typically contain both a numeric `Well`
column (`1, 2, 3 …`) and an alphanumeric `Well Position` column
(`A1, A2 …`). The detector picks `Well Position` automatically because it
recognises the A1-format well IDs.

### Supported file formats

| Extension      | Parser                              |
|----------------|-------------------------------------|
| `.xlsx`        | `openpyxl`                          |
| `.xls`         | `xlrd` (install `qpcr-analyzer[xls]`) |
| `.csv`         | comma-separated                     |
| `.tsv`, `.txt` | tab-separated                       |

---

## Methods

### ΔCt

```
dCt(sample, target) = mean_Cq(target, sample) − mean_Cq(HK, sample)
Expr_vs_HK          = 2^(−dCt)
```

Quantifies expression relative to the housekeeping gene only. No reference
biological group is involved. Useful for comparing individual sample-level
HK-normalised expression or for datasets without a clear control group.

### ΔΔCt (batch-aware)

```
dCt(sample, target)    = mean_Cq(target, sample) − mean_Cq(HK, sample)

Ref_dCt(batch, target) = mean of dCt over all samples in the
                         reference group that belong to this batch

ddCt(sample, target)   = dCt(sample, target) − Ref_dCt(batch(sample), target)

Relative_Expr          = 2^(−ddCt)
```

Each batch is normalised independently so that the **reference group's mean
ΔΔCt = 0 within that batch** (the exact mathematical invariant). After
normalisation, all batches are concatenated. This eliminates run-to-run
variation while preserving biological differences.

> Note: `mean(2^(−ΔΔCt))` of the reference group is **not** exactly 1
> unless reference samples have identical ΔCt — exponentials don't commute
> with averaging. The anchored quantity is the additive ΔΔCt, not the
> multiplicative relative expression.

**Requirement:** every batch must contain at least one sample from the
reference group for every measured target gene.

---

## Output workbook

The downloaded `.xlsx` file contains:

| Sheet              | Contents                                                      |
|--------------------|---------------------------------------------------------------|
| `raw_data`         | Standardised well-level data, with `Excluded` flag.           |
| `dCt_{HK}`         | Full ΔCt table (one per housekeeping gene).                   |
| `ddCt_{HK}`        | Full ΔΔCt table including `Batch` and `Reference_Group`.      |
| `prism_dCt_{HK}`   | Wide-format ΔCt table for GraphPad Prism.                     |
| `prism_ddCt_{HK}`  | Wide-format Relative_Expr table for GraphPad Prism.           |

### Prism sheet format

```
Target         Control    Treatment    KO
GeneX
               1.05       2.31         0.45
               0.93       2.15         0.52
               1.02       2.48         0.41

GeneY
               0.88       1.72         …
```

Each block is a separate target gene. Paste directly into a Prism
**Grouped** table (rows = replicates, column sets = groups).

---

## Programmatic API

The `qpcr_analyzer.core` package can be used as a library, independent of
the UI:

```python
import pandas as pd
from qpcr_analyzer.core import (
    read_table,
    detect_columns,
    apply_mapping,
    validate_sample_groups,
    mark_outliers,
    compute_delta_ct,
    compute_delta_delta_ct,
    results_to_xlsx_bytes,
)

# 1. Load data
df = read_table("my_experiment.xlsx", "my_experiment.xlsx")

# 2. Detect and apply column mapping
mapping = detect_columns(df)
assert mapping.validate() == []                # [] if all required columns found
std = apply_mapping(df, mapping)

# 3. Validate sample ↔ group consistency
errors = validate_sample_groups(std)
if errors:
    raise ValueError("\n".join(errors))

# 4. Mark outlier wells (tolerance = 1 Cq cycle)
flagged = mark_outliers(std, tolerance=1.0)
std["Excluded"] = flagged["Outlier"]

# 5. ΔCt
dct = compute_delta_ct(std, ref_genes=["GAPDH"])

# 6. Batch-aware ΔΔCt
batches = {"s1": "run1", "s2": "run1", "s3": "run2", "s4": "run2"}
ddct = compute_delta_delta_ct(
    std,
    ref_genes=["GAPDH"],
    reference_group="Control",
    sample_batches=batches,        # omit for single-batch data
)

# 7. Export
xlsx = results_to_xlsx_bytes(std, dct, ddct)
with open("results.xlsx", "wb") as f:
    f.write(xlsx)
```

---

## Project layout

```
qpcr-analyzer/
├── pyproject.toml                      build config & dependencies
├── README.md
├── LICENSE
├── src/qpcr_analyzer/
│   ├── __init__.py                     __version__
│   ├── __main__.py                     console-script entry point
│   ├── app/
│   │   ├── __init__.py
│   │   └── main.py                     NiceGUI 5-step stepper UI + start()
│   └── core/                           pure-Python, no UI deps
│       ├── __init__.py                 re-exports the public API
│       ├── io.py                       file readers
│       ├── columns.py                  column-role detection & mapping
│       ├── outliers.py                 replicate outlier flagging
│       ├── quant.py                    ΔCt and batch-aware ΔΔCt
│       └── export.py                   Excel writer (raw + ΔCt + ΔΔCt + Prism)
└── tests/
    ├── test_columns.py
    ├── test_outliers.py
    └── test_quant.py
```

The split between `core/` and `app/` is deliberate:

* `core/` is fully unit-testable without a browser.
* `app/main.py` is a thin presentation layer — every scientific decision
  happens in `core/` and is therefore covered by the test suite.

---

## Testing

```bash
pip install -e ".[dev]"
pytest                          # all unit tests
pytest --cov=qpcr_analyzer       # with coverage
ruff check .                    # lint
```

The test suite covers column detection (including the Applied Biosystems
"Well vs Well Position" tiebreak), outlier flagging across replicate
counts, ΔCt arithmetic, batch-aware ΔΔCt anchoring, and reference-group
validation. UI code is intentionally not tested — it is a thin shell.

---

## License

MIT © 2026 Jielin Yang
