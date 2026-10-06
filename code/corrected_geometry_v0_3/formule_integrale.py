"""Traceable C estimate: finite large-a integral, then a labelled model remainder.

Never interpret mc_se as a complete uncertainty or this estimate as a proof.
Independent draws are used at quadrature nodes. The targeted window may omit
contributions; neither that omission nor asymptotic tails have certified bounds.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import scipy
from integral_model import echantillon, G_cible

G0 = 4 * math.sqrt(2) / 15
SIEGEL = 6 / math.pi ** 2

def log_nodes(lo, hi, intervals, order):
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError('bornes logarithmiques invalides')
    if (not isinstance(intervals, int) or not isinstance(order, int)
            or isinstance(intervals, bool) or isinstance(order, bool)
            or intervals < 1 or order < 1):
        raise ValueError('nombres entiers de panneaux/noeuds positifs requis')
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(lo, hi, intervals + 1)
    for u0, u1 in zip(edges[:-1], edges[1:]):
        for xx, ww in zip(x, w):
            u = (u0 + u1) / 2 + (u1 - u0) / 2 * xx
            a = math.exp(u)
            yield a, (u1 - u0) / 2 * float(ww) * a * a

def deterministic_tail(G, A, Amax, intervals=4, order=6):
    """Integral of a*G(a) over [A,Amax], without replacing a^3 G by a mean."""
    return sum(weight * G(a) for a, weight in log_nodes(math.log(A), math.log(Amax), intervals, order))

def estimate(*, seed=42012, n=10000, target_n=10000, intervals=12, order=4,
             tail_intervals=2, tail_order=6, umin=-12.0, A=6.0, Amax=16.0, window=5.0):
    if not (math.isfinite(A) and math.isfinite(Amax) and 0 < A < Amax
            and umin < math.log(A) and 0 < window and 2 * window < A ** 3):
        raise ValueError('configuration de queues invalide')
    if not (isinstance(n, int) and isinstance(target_n, int) and n >= 2 and target_n >= 2):
        raise ValueError('tailles Monte Carlo entieres >=2 requises')
    rng = np.random.default_rng(seed)
    middle = tail = variance = Nmiddle = Nvariance = 0.0
    rows = []
    for a, weight in log_nodes(umin, math.log(A), intervals, order):
        c, ar = echantillon(a, n, rng)
        mean, se = float(c.mean()), float(c.std(ddof=1) / math.sqrt(n))
        middle += weight * mean
        variance += (weight * se) ** 2
        Nmiddle += weight * float(ar.mean())
        Nvariance += (weight * float(ar.std(ddof=1)) / math.sqrt(n)) ** 2
        rows.append(('middle_uniform', a, weight, mean, se, n))
    for a, weight in log_nodes(math.log(A), math.log(Amax), tail_intervals, tail_order):
        mean, se = G_cible(a, target_n, rng, window, window)
        tail += weight * mean
        variance += (weight * se) ** 2
        rows.append(('large_finite_target_window', a, weight, float(mean), float(se), target_n))
    small_model = 2 * G0 * math.exp(umin / 2)
    large_model = 1 / (9 * Amax)
    result = {
        'status': 'NUMERICAL_ESTIMATE_NOT_CERTIFIED',
        'C_estimate': SIEGEL * (middle + tail + small_model + large_model),
        'C_mc_se_only': SIEGEL * math.sqrt(variance),
        'components_C': {'middle': SIEGEL * middle, 'large_finite_window': SIEGEL * tail,
                         'small_asymptotic_model': SIEGEL * small_model,
                         'large_asymptotic_model': SIEGEL * large_model},
        'N_truncated_at_A': SIEGEL * Nmiddle,
        'N_truncated_mc_se_only': SIEGEL * math.sqrt(Nvariance),
        'unbounded_errors': ['small-tail asymptotic remainder', 'large-tail asymptotic remainder',
                             'target-window coverage', 'quadrature discretization',
                             'floating-point membership'],
        'config': dict(seed=seed, n=n, target_n=target_n, intervals=intervals, order=order,
                       tail_intervals=tail_intervals, tail_order=tail_order,
                       umin=umin, A=A, Amax=Amax, window=window),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'platform': platform.platform()},
        'code_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted(Path(__file__).parent.glob('*.py'))},
    }
    return result, rows

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile', choices=['smoke', 'validation', 'reference'], default='smoke')
    p.add_argument('--seed', type=int, default=42012)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    p.error("Recalcul de C suspendu dans v0.3; audit geometrique uniquement")
    settings = {
        'smoke': {},
        'validation': dict(n=200000, target_n=50000, intervals=29, order=6, tail_intervals=3, tail_order=6),
        'reference': dict(n=1000000, target_n=200000, intervals=29, order=6, tail_intervals=3, tail_order=6),
    }[args.profile]
    # Append files only to a new output directory; never overwrite an earlier estimate.
    args.out.mkdir(parents=True, exist_ok=False)
    result, rows = estimate(seed=args.seed, **settings)
    (args.out / 'result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    with (args.out / 'nodes.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['channel', 'a', 'weight_a_da', 'G_mean', 'G_mc_se', 'n'])
        writer.writerows(rows)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
