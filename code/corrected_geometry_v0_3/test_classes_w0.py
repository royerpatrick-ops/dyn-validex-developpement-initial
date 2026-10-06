"""Disk edge classes: corrected g0 and uncertainty including absent classes.

This is a finite-R comparison to a model asymptotic, not a proof or PASS criterion.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull
from geometry_columns import column_indices, integer_hull_vertices

K = 4 * math.sqrt(2) / 15

def primitive_counts(max_r2):
    bound = math.isqrt(max_r2)
    counts = {}
    for x in range(-bound, bound + 1):
        for y in range(-bound, bound + 1):
            r2 = x*x + y*y
            if 0 < r2 <= max_r2 and math.gcd(x, y) == 1:
                counts[r2] = counts.get(r2, 0) + 1
    return counts

def classes(R, N, rng, max_r2=100):
    if not math.isfinite(R) or R < 2 or not isinstance(N, int) or N < 2:
        raise ValueError('R >= 2 et N entier >=2 requis')
    if not isinstance(max_r2, int) or max_r2 < 1:
        raise ValueError('max_r2 entier positif requis')
    counts = primitive_counts(max_r2)
    samples = {k: np.zeros(N) for k in counts}
    for i in range(N):
        th = rng.uniform(0, math.pi/2)
        u, v = rng.uniform(0, 1, 2)
        indices = column_indices('ellipse', 1.0, float(R), th, u, v)
        w = integer_hull_vertices(indices).astype(float) - [u, v]
        e = np.roll(w, -1, axis=0) - w
        angles = np.arctan2(w[:,1], w[:,0])
        arcs = np.mod(np.roll(angles, -1) - angles, 2*math.pi)
        c = R * arcs - np.hypot(e[:,0], e[:,1])
        ei = np.rint(e).astype(np.int64)
        if np.max(np.abs(e - ei)) > 1e-7:
            raise ArithmeticError('aretes non entieres')
        g = np.gcd(np.abs(ei[:,0]), np.abs(ei[:,1]))
        if np.any(g == 0):
            raise ArithmeticError('arete nulle')
        r2 = (ei[:,0]//g)**2 + (ei[:,1]//g)**2
        for k in counts:
            samples[k][i] = c[r2 == k].sum() * R**(1/3) / (2*math.pi)
    return [{'r2': k, 'n_primitive': counts[k], 'mean': float(v.mean()),
             'mc_se': float(v.std(ddof=1)/math.sqrt(N)),
             'model_prediction': counts[k]*K*R**(-1/6)*k**(-.75)/(2*math.pi)}
            for k, v in sorted(samples.items())]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--R', type=float, default=10000)
    p.add_argument('--samples', type=int, default=16)
    p.add_argument('--seed', type=int, default=42013)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    p.error("Tests Monte Carlo par classes suspendus dans v0.3; audit geometrique uniquement")
    result = {'status':'FINITE_R_MODEL_COMPARISON', 'R':args.R, 'samples':args.samples,
              'seed':args.seed, 'K':K,
              'classes':classes(args.R, args.samples, np.random.default_rng(args.seed))}
    with args.out.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
