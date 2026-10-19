"""Tests of the PmagPy field model and sun compass functions used by mk_sam_file.py.

Run with:  python -m pytest tests
"""

import csv
import os
import shutil
import subprocess
import sys

import pytest
from pmagpy import ipmag, pmag

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO_DIR, 'tests', 'data')


def test_pmagpy_uses_igrf14():
    """The installed PmagPy has IGRF-14 (pmagpy >= 4.3.13).

    At the N7 Graafwater sections in May 2025 IGRF-14 gives a declination of -24.59
    degrees, consistent with the sun-minus-magnetic declinations of the GW1-GW4 sun
    compass orientations (site medians -24.2 to -25.6). Older PmagPy versions
    (IGRF-13) give a different value.
    """
    dec, inc, intensity = ipmag.igrf([2025.41, 0.0, -31.975637, 18.741023])
    assert dec - 360 == pytest.approx(-24.59, abs=0.01)
    assert inc == pytest.approx(-65.18, abs=0.01)


def test_sun_compass_declination():
    """pmag.dosundec reproduces the GW3-1 sun compass core strike (147.7 degrees)."""
    sundata = {'date': '2025:05:31:15:52', 'delta_u': 2, 'lat': '-31.974709',
               'lon': '18.740605', 'shadow_angle': 374.0}
    assert round(pmag.dosundec(sundata), 1) == pytest.approx(147.7)


def test_dates_after_igrf_model_raise(tmp_path):
    """A sampling date after 2030 stops mk_sam_file.py instead of using a zero field."""
    csv_path = tmp_path / 'GW3-.csv'
    with open(os.path.join(DATA_DIR, 'GW3-', 'input', 'GW3-.csv'), newline='') as f:
        rows = list(csv.reader(f))
    header_index = next(i for i, row in enumerate(rows) if row and row[0] == 'sample_name')
    year_column = rows[header_index].index('year')
    for row in rows[header_index + 1:]:
        if row and row[0]:
            row[year_column] = '2031'
    with open(csv_path, 'w', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(rows)
    result = subprocess.run([sys.executable, 'mk_sam_file.py', str(csv_path), str(tmp_path)],
                            cwd=REPO_DIR, capture_output=True, text=True)
    assert result.returncode != 0
    assert 'after the end of the IGRF-14 model' in result.stderr
