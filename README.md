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
6. [Outputs](#outputs)
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

### Install from PyPI (recommended)

Once a virtual environment is active (see per-OS setup below):

```bash
pip install qpcr-analyzer                    # core install
pip install "qpcr-analyzer[xls]"             # + legacy .xls support
qpcr-analyzer                                # launch the app
```

#### Per-OS setup (Python + virtual environment)

<details>
<summary><b>Windows</b></summary>

PowerShell, from any folder you like:

```powershell
# 1. Install Python 3.10+ from https://www.python.org/downloads/windows/
#    (tick "Add python.exe to PATH" in the installer)

# 2. Create and activate a virtual environment
python -m venv qpcr-venv
qpcr-venv\Scripts\Activate.ps1

# 3. Install qPCR Analyzer from PyPI
pip install qpcr-analyzer

# 4. Run it
qpcr-analyzer
```

If PowerShell blocks the activation script, run once:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

</details>

<details>
<summary><b>macOS</b></summary>

```bash
# 1. Install Python (Homebrew is easiest)
brew install python@3.12

# 2. Create and activate a virtual environment
python3 -m venv qpcr-venv
source qpcr-venv/bin/activate

# 3. Install qPCR Analyzer from PyPI
pip install qpcr-analyzer

# 4. Run it
qpcr-analyzer
```

On Apple Silicon, all required wheels (`pandas`, `numpy`, `openpyxl`,
`plotly`, `rapidfuzz`, `nicegui`) are published as native arm64, so
installation does not need to compile from source.

</details>

<details>
<summary><b>Linux</b></summary>

```bash
# 1. Make sure python3-venv is available (Debian/Ubuntu)
sudo apt install python3 python3-venv python3-pip      # Debian/Ubuntu
# sudo dnf install python3 python3-pip                  # Fedora/RHEL
# sudo pacman -S python python-pip                      # Arch

# 2. Create and activate a virtual environment
python3 -m venv qpcr-venv
source qpcr-venv/bin/activate

# 3. Install qPCR Analyzer from PyPI
pip install qpcr-analyzer

# 4. Run it
qpcr-analyzer
```

</details>

#### Optional extras

| Extra  | Purpose | Install command |
|--------|---------|-----------------|
| `xls`  | Read legacy `.xls` files (Excel 97-2003) | `pip install "qpcr-analyzer[xls]"` |
| `dev`  | Test runner + linter for contributors | `pip install -e ".[dev]"` (after cloning) |

#### Installing a specific version

```bash
pip install qpcr-analyzer==2.1.0
```

#### Upgrading

```bash
pip install --upgrade qpcr-analyzer
```

#### Uninstalling

```bash
pip uninstall qpcr-analyzer
```

### Install from source (GitHub)

For unreleased changes from `main`, or to develop against the codebase:

```bash
# Latest commit on main
pip install "git+https://github.com/j-y26/py_qpcr_analyzer.git"

# A specific tag
pip install "git+https://github.com/j-y26/py_qpcr_analyzer.git@v2.1.0"

# With the xls extra
pip install "qpcr-analyzer[xls] @ git+https://github.com/j-y26/py_qpcr_analyzer.git"
```

---

## Running the application

### Test environment (development/trial)

Use this when you want to try the app, run the test suite, or modify the
source code.

```bash
git clone https://github.com/j-y26/py_qpcr_analyzer.git
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
/opt/qpcr-analyzer/bin/pip install "qpcr-analyzer[xls]"

# 2. Run the server bound to all interfaces on the standard port
QPCR_HOST=0.0.0.0 QPCR_PORT=8090 /opt/qpcr-analyzer/bin/qpcr-analyzer
```

Then point colleagues at `http://<server-hostname>:8090`.

#### Configuration (environment variables)

| Variable    | Default     | Description                                          |
|-------------|-------------|------------------------------------------------------|
| `QPCR_HOST` | `127.0.0.1` | Bind address. Use `0.0.0.0` to expose on the LAN.    |
| `QPCR_PORT` | `8090`      | TCP port. |

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
  unless the user clicks one of the **Download** buttons (xlsx, CSV
  bundle, or per-figure PNG).

---

## Usage walkthrough

The left pane is a seven-step stepper; each step's *Continue* button
unlocks once the prior step is complete. The right pane shows live
results across tabs (*Summary*, *Data preview*, *Excluded blocks*,
*Housekeeping & exclusions*, *ΔCt results*, *ΔΔCt results*, *Downloads*).

| Step | What you do |
|------|-------------|
| **1. Upload data**            | Drop in an `.xlsx`, `.xls`, `.csv`, `.tsv`, or `.txt` file (one row per well). |
| **2. Column mapping**         | Confirm or adjust auto-detected role assignments (Well, Target, Sample, Cq, plus optional Group and Batch). |
| **3. Groups & batches**       | Review or assign each sample's biological group. Optionally toggle multi-batch mode and assign samples to **batches** if they come from multiple independent runs. |
| **4. Outliers**               | Adjust replicate tolerance and exclude outlier wells. The app keeps the tightest contiguous run of replicates within tolerance and flags the rest; NaN Cq values are always flagged. |
| **5. Housekeeping gene(s)**   | Pick one or more housekeeping genes, review which samples lack a usable HK Cq (skipped per HK), and click **Apply** to lock the selection. |
| **6. Run ΔCt**                | Compute ΔCt per housekeeping gene. No reference group needed — bar plots and downloads populate immediately. |
| **7. Run ΔΔCt** *(optional)*  | Pick the **reference group** (ΔΔCt anchor) and run batch-aware ΔΔCt. Each batch is anchored independently; results are merged for plotting. |

### Column auto-detection

The detector scores every column name against a synonym list using exact,
substring, and fuzzy matching (`rapidfuzz`). Columns scoring below 0.85
are left unassigned and shown as "unmatched" in the UI. **Group** and
**Batch** are optional — leaving them as `(none)` is fine, you assign them
in step 3.

| Role   | Required | Example column names                              |
|--------|----------|---------------------------------------------------|
| Well   | yes      | Well, Well Position, Location                     |
| Target | yes      | Target, Gene, Assay, Detector                     |
| Sample | yes      | Sample, Sample ID, Sample Name                    |
| Cq     | yes      | Cq, Ct, Cq Value, Ct Mean                         |
| Group  | no       | Group, Condition, Treatment, Biological Set Name  |
| Batch  | no       | Batch, Run, Plate, Experiment                     |

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

## Outputs

The **Downloads** tab offers three artefacts. All three reflect whatever
exclusions, housekeeping selections, and reference group are currently
configured — re-run ΔCt / ΔΔCt after any change to refresh them.

### Excel workbook (`.xlsx`)

| Sheet                  | Contents                                                                                                          |
|------------------------|-------------------------------------------------------------------------------------------------------------------|
| `raw_data`             | Standardised well-level data, sorted Target → Sample → Well, with `Excluded` flag.                                |
| `dCt_{HK}`             | Long-format ΔCt table per housekeeping gene. Columns: `Target, Group, Sample, Mean_Cq, Mean_Cq_{HK}, dCt, Expr_vs_HK, Reference_Gene`. |
| `formatted_dCt_{HK}`   | Wide grouped table of **`Expr_vs_HK` (= 2^−ΔCt)** with sample names per row — see layout below.                   |
| `ddCt_{HK}`            | Long-format ΔΔCt table per housekeeping gene. Columns: `Target, Group, Sample, Batch, Mean_Cq, Mean_Cq_{HK}, dCt, Ref_dCt, ddCt, Relative_Expr, Is_Reference_Group, Reference_Gene, Reference_Group`. |
| `formatted_ddCt_{HK}`  | Wide grouped table of **`Relative_Expr` (= 2^−ΔΔCt)** with sample names per row.                                  |

The `Batch` column is included on `ddCt_*` sheets only when multi-batch
mode is enabled in step 3; single-batch exports omit it.

### CSV bundle (`.zip`)

The same logical sheets as the Excel workbook, packed as one CSV per
file inside a single zip archive (`raw_data.csv`, `dCt_*.csv`,
`ddCt_*.csv`, `formatted_*.csv`). Use this when downstream tools prefer
plain text.

### Figures (`.png`)

Each individual ΔCt / ΔΔCt bar plot has a built-in PNG export in its
toolbar (camera icon). The **Downloads** tab additionally exposes one
combined figure per housekeeping gene — every target packed into a
single PNG — for ΔCt and (if run) ΔΔCt.

### Formatted sheet layout

```
Target  | Sample (Control) | Relative Expression (Control) | Sample (Treatment) | Relative Expression (Treatment) | …
GeneX   |                  |                               |                    |                                 |
        | s1               | 1.05                          | s4                 | 2.31                            |
        | s2               | 0.93                          | s5                 | 2.15                            |
        | s3               | 1.02                          | s6                 | 2.48                            |
        |                  |                               |                    |                                 |
GeneY   |                  |                               |                    |                                 |
        | s1               | 0.88                          | s4                 | 1.72                            |
        | …                | …                             | …                  | …                               |
```

One block per target gene, separated by a blank row. Each biological
group contributes two adjacent columns — the sample name and its value —
so row provenance is preserved. Groups with fewer samples are padded
with blanks. Every value is parity-checked against the long-format
source at build time; a mismatch raises `ValueError` rather than
silently writing inconsistent data.

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
    results_to_csv_zip_bytes,
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

# 7. Export — Excel workbook (raw + per-HK ΔCt/ΔΔCt + formatted_*)
xlsx = results_to_xlsx_bytes(std, dct, ddct)
with open("results.xlsx", "wb") as f:
    f.write(xlsx)

# Or the same logical sheets as a zip of CSVs
csv_zip = results_to_csv_zip_bytes(std, dct, ddct)
with open("results.zip", "wb") as f:
    f.write(csv_zip)
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
│   │   └── main.py                     NiceGUI 7-step stepper UI + start()
│   └── core/                           pure-Python, no UI deps
│       ├── __init__.py                 re-exports the public API
│       ├── io.py                       file readers
│       ├── columns.py                  column-role detection & mapping
│       ├── outliers.py                 replicate outlier flagging
│       ├── quant.py                    ΔCt and batch-aware ΔΔCt
│       ├── summary.py                  well/sample/target ordering helpers
│       └── export.py                   xlsx + CSV-zip writers (raw + ΔCt + ΔΔCt + formatted)
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
