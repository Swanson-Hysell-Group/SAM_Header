from __future__ import annotations

import os
from pathlib import Path

import PyInstaller.__main__


ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / 'GUI' / 'SAMHeaderBuilder.spec'


def main() -> None:
    os.chdir(ROOT)
    PyInstaller.__main__.run(
        [
            '--noconfirm',
            '--clean',
            '--distpath',
            str(ROOT / 'GUI' / 'dist'),
            '--workpath',
            str(ROOT / 'GUI' / 'build'),
            str(SPEC),
        ]
    )


if __name__ == '__main__':
    main()