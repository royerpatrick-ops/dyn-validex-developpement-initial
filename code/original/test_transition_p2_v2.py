import numpy as np

import transition_p2_v2 as mod


def test_fit_alpha_gamma_ignores_nonpositive_residuals():
    R = np.array([100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0])
    resid = np.array([
        5.0e-3,
        4.0e-3,
        3.0e-3,
        2.0e-3,
        1.0e-3,
        -1.0e-3,
        -2.0e-3,
    ])

    alpha, gamma, r2 = mod.fit_alpha_gamma(R, resid, R_min=100, R_max=400)

    assert np.isfinite(alpha)
    assert np.isfinite(gamma)
    assert np.isfinite(r2)
    assert alpha > 0
    assert gamma > 0


def test_fit_alpha_gamma_uses_absolute_magnitude_when_near_transition():
    R = np.array([100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 350.0, 400.0])
    resid = np.array([
        -5.0e-3,
        -4.0e-3,
        -3.0e-3,
        -2.0e-3,
        -1.5e-3,
        -1.0e-3,
        -7.5e-4,
        -5.0e-4,
        -3.5e-4,
        -2.0e-4,
        -1.0e-4,
    ])

    alpha, gamma, r2 = mod.fit_alpha_gamma(R, resid, R_min=100, R_max=400)

    assert np.isfinite(alpha)
    assert np.isfinite(gamma)
    assert np.isfinite(r2)
    assert alpha > 0
    assert gamma > 0
