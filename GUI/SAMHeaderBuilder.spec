# -*- mode: python ; coding: utf-8 -*-

import os
from pathlib import Path

ROOT = Path(os.getcwd()).resolve()

a = Analysis(
    [str(ROOT / 'GUI' / 'app.py')],
    pathex=[str(ROOT)],
    binaries=[],
    datas=[
        (str(ROOT / 'GUI' / 'assets' / 'sam_header_icon.svg'), 'GUI/assets'),
        (str(ROOT / 'GUI' / 'assets' / 'sam_header_icon.ico'), 'GUI/assets'),
    ],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='SAMHeaderBuilder',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ROOT / 'GUI' / 'assets' / 'sam_header_icon.ico'),
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='SAMHeaderBuilder',
)
