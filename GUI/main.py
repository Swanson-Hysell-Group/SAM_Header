from __future__ import annotations

import os
import sys
from pathlib import Path


def _venv_python_candidate(repo_root: Path) -> Path | None:
    candidates = [
        repo_root / 'venv' / 'Scripts' / 'python.exe',
        repo_root / '.venv' / 'Scripts' / 'python.exe',
        repo_root / 'venv' / 'bin' / 'python',
        repo_root / '.venv' / 'bin' / 'python',
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _maybe_reexec_with_venv() -> None:
    # Prevent recursion if the selected interpreter cannot import dependencies.
    if os.environ.get('SAM_HEADER_GUI_REEXEC') == '1':
        return

    repo_root = Path(__file__).resolve().parents[1]
    venv_python = _venv_python_candidate(repo_root)
    if venv_python is None:
        return

    current_python = Path(sys.executable).resolve()
    if current_python == venv_python.resolve():
        return

    os.environ['SAM_HEADER_GUI_REEXEC'] = '1'
    os.execv(str(venv_python), [str(venv_python), str(Path(__file__).resolve()), *sys.argv[1:]])


_maybe_reexec_with_venv()

from app import main


if __name__ == '__main__':
    raise SystemExit(main())