"""Numerical curvature integrals with analytic finiteness checks.

For |x|^p+|y|^p=1, p>1, near an axis: kappa^a ds ~ x^(a(p-2)) dx.
Thus I_a is finite iff 1+a(p-2)>0. Numerical error estimates are not bounds.
"""
import math
import warnings
from scipy.integrate import quad, IntegrationWarning

class DivergentInvariant(ValueError):
    pass

def parameters(kind, prm, a=0.0):
    if kind not in ('super', 'ellipse'):
        raise ValueError('forme inconnue')
    if not math.isfinite(prm) or not math.isfinite(a) or a < 0:
        raise ValueError('parametre fini et exposant fini non negatif requis')
    if kind == 'ellipse' and prm <= 0:
        raise ValueError('rapport des axes positif requis')
    if kind == 'super':
        if prm <= 1:
            raise ValueError('superellipse: p > 1 requis')
        margin = 1 + a * (prm - 2)
        # Conservative rejection at the floating representation of the critical exponent.
        if margin <= 16 * math.ulp(1.0):
            raise DivergentInvariant(f'I_{a:g} non fini ou seuil numeriquement critique: p={prm:g}')

def checked_quad(f, lo, hi):
    with warnings.catch_warnings():
        warnings.simplefilter('error', IntegrationWarning)
        value, error = quad(f, lo, hi, limit=2000, epsabs=1e-12, epsrel=2e-12)
    if not (math.isfinite(value) and math.isfinite(error) and value > 0 and error >= 0):
        raise ArithmeticError('quadrature non finie ou non positive')
    if error > max(1e-10, abs(value) * 1e-9):
        raise ArithmeticError('quadrature insuffisamment resolue')
    return value, error

def super_result(p, a):
    parameters('super', p, a)
    margin = 1 + a * (p - 2)
    m = max(6.0, 2.0 / margin)
    xm = 2 ** (-1 / p)
    # x=xm*t^m; collect powers before evaluation to avoid singular factors.
    factor = (p - 1) ** a * xm ** margin * m
    def f(t):
        x = xm * t ** m
        w = 1 - x ** p
        slope2 = x ** (2 * (p - 1)) * w ** (2 * (1 / p - 1))
        return (factor * t ** (m * margin - 1) * w ** (a * (1 / p - 2))
                * (1 + slope2) ** ((1 - 3 * a) / 2))
    v, e = checked_quad(f, 0.0, 1.0)
    return {'value': 8 * v, 'quad_error_estimate': 8 * e, 'certified': False}

def ellipse_result(k, a):
    parameters('ellipse', k, a)
    b = 1 / k
    def f(t):
        speed = math.hypot(math.sin(t), b * math.cos(t))
        return (b / speed ** 3) ** a * speed
    v, e = checked_quad(f, 0.0, math.pi / 2)
    return {'value': 4 * v, 'quad_error_estimate': 4 * e, 'certified': False}

def I_super(p, a):
    return super_result(p, a)['value']

def I_ellipse(k, a):
    return ellipse_result(k, a)['value']

def invariants(kind, prm):
    parameters(kind, prm)
    f = I_super if kind == 'super' else I_ellipse
    return f(prm, 0), f(prm, 1 / 3), f(prm, 4 / 3)
