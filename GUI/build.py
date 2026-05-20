from __future__ import annotations

from pathlib import Path

import PyInstaller.__main__


ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / 'GUI' / 'app.py'


def main() -> None:
    PyInstaller.__main__.run(
        [
            '--noconfirm',
            '--clean',
            '--windowed',
            '--name',
            'SAMHeaderBuilder',
            '--distpath',
            str(ROOT / 'GUI' / 'dist'),
            '--workpath',
            str(ROOT / 'GUI' / 'build'),
            '--specpath',
            str(ROOT / 'GUI'),
            str(APP),
        ]
    )


if __name__ == '__main__':
    main()