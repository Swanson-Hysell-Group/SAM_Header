"""Tests of the IGRF-14 field calculation used to correct magnetic compass orientations.

Reference values were calculated with pmagpy.pmag.doigrf (PmagPy commit 3752bf41),
which uses the same IGRF-14 coefficient file and the same synthesis routine (magsyn).

Run with:  python -m pytest tests
"""

import hashlib
import os
import sys

import numpy as np
import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_DIR)

import coefficients as cf  # noqa: E402
import mk_sam_utilities as sam  # noqa: E402

# (lon, lat, alt_km, date, x_nT, y_nT, z_nT) from pmagpy.pmag.doigrf
PMAGPY_REFERENCE = [
    (18.741023, -31.975637, 0.0, 2025.41, 9688.4431, -4432.9722, -23034.5589),
    (18.966727, -32.250497, 0.0, 2026.8, 9699.5096, -4583.8934, -22900.5195),
    (-93.265, 44.9778, 0.25, 1987.3, 16997.153, 1011.3644, 55913.5424),
    (151.2, -33.9, 0.0, 1950.0, 25113.952, 4753.7055, -52390.5736),
    (-70.0, 80.0, 1.0, 2029.5, 3063.4378, -2309.7353, 56075.0224),
    (0.0, 0.0, 0.0, 1912.6, 28080.8337, -8197.1226, -7369.0661),
]

# SHA-256 of https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt (downloaded 10/19/2026)
IGRF14_SHA256 = '8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c'


def test_coefficient_file_is_official_igrf14():
    """The bundled coefficient file is unmodified from the IAGA distribution."""
    path = os.path.join(REPO_DIR, 'field_models', 'igrf14coeffs.txt')
    with open(path, 'rb') as f:
        assert hashlib.sha256(f.read()).hexdigest() == IGRF14_SHA256


def test_igrf14_epochs_and_coefficients():
    """Epochs run 1900-2030 every 5 years, and spot values match the IGRF-14 table."""
    models, coeffs = cf.get_igrf14()
    assert models == list(range(1900, 2031, 5))
    assert coeffs.shape == (27, 195)
    assert not np.isnan(coeffs).any()
    g10, g11, h11 = coeffs[models.index(2025)][:3]
    assert (g10, g11, h11) == (-29350.0, -1410.3, 4545.5)
    # 2030 epoch = 2025 main field + 5 years of 2025-2030 secular variation (g10 SV = 12.6 nT/yr)
    assert coeffs[models.index(2030)][0] == pytest.approx(-29350.0 + 5 * 12.6)


@pytest.mark.parametrize('lon, lat, alt, date, x, y, z', PMAGPY_REFERENCE)
def test_doigrf_matches_pmagpy(lon, lat, alt, date, x, y, z):
    """Field components agree with pmagpy.pmag.doigrf to 0.001 nT."""
    result = sam.doigrf(lon, lat, alt, date)
    assert np.allclose(result[:3], [x, y, z], atol=1e-3)


def test_field_at_epoch_uses_epoch_coefficients():
    """At an epoch the field depends only on that epoch's coefficients."""
    models, coeffs = cf.get_igrf14()
    gh = coeffs[models.index(2020)]
    expected = sam.magsyn(gh, np.zeros(len(gh)), 2020, 2020.0, 1, 0.0, 90 + 31.975637, 18.741023)
    assert np.allclose(sam.doigrf(18.741023, -31.975637, 0.0, 2020.0), expected)


def test_interpolation_is_linear_in_coefficients():
    """Midway between epochs the field is synthesized from the mean coefficients."""
    models, coeffs = cf.get_igrf14()
    gh_mid = (coeffs[models.index(2015)] + coeffs[models.index(2020)]) / 2
    expected = sam.magsyn(gh_mid, np.zeros(len(gh_mid)), 2017.5, 2017.5, 1, 0.0,
                          90 + 31.975637, 18.741023)
    assert np.allclose(sam.doigrf(18.741023, -31.975637, 0.0, 2017.5), expected, atol=1e-6)


def test_declination_at_graafwater_sites_2025():
    """IGRF-14 declination at the N7 Graafwater sections in May 2025 is about -24.6 degrees.

    This agrees with the sun-minus-magnetic declinations of the GW1-GW4 sun compass
    orientations (site medians -24.2 to -25.6 degrees).
    """
    dec, inc, intensity = sam.igrf([2025.41, 0.0, -31.975637, 18.741023])
    assert dec - 360 == pytest.approx(-24.59, abs=0.01)


def test_dates_after_2030_raise():
    with pytest.raises(ValueError):
        sam.doigrf(18.7, -32.0, 0.0, 2030.5)
