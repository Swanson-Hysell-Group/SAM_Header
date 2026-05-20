# SAM Header GUI draft

This folder contains a PySide6 desktop wrapper around the repository's CSV-to-SAM conversion logic.

## What it does

- loads the existing CSV template workflow
- keeps Windows-style CRLF output for generated files on Windows and macOS
- lets users choose a core or block workflow
- exposes separate declination toggles for core strike, block strike, and bedding strike
- packages into a standalone desktop app with PyInstaller

## Local setup

From the repository root:

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r GUI/requirements.txt
```

## Run the draft app

```bash
python GUI/app.py
```

## Build a desktop app

```bash
python GUI/build.py
```

PyInstaller outputs are written to `GUI/dist` and intermediate files to `GUI/build`.