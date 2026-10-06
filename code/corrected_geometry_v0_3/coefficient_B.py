"""Internal B-model check; no claim about the hull's second-order asymptotic.

Integrate rho on [1/T,infinity) by finite primitive-vector sums, row by row.
The cutoff convergence and Mellin identification remain model calculations.
"""
import argparse
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np

def cutoff_integral(T):
    if not isinstance(T, int) or T < 2:
        raise ValueError('T entier >=2 requis')
    p = np.arange(-T, T+1, dtype=np.int64)
    sums = np.zeros(3)
    for q in range(T+1):
        pp = p if q else p[p > 0]
        r2 = pp*pp + q*q
        mask = (r2 < T*T) & (np.gcd(pp, q) == 1)
        r = np.sqrt(r2[mask].astype(float))
        sums += 2 * np.array([np.sum(r**-1.5), np.sum(1/r), np.sum(r)])
    lam = 1/T
    return math.pi/8 * ((8/5)*sums[0] - 2*math.sqrt(lam)*sums[1]
                        + (2/5)*lam**2.5*sums[2]) - 2/math.sqrt(lam)

def model_coefficients():
    with mp.workdps(40):
        Z = 4*mp.zeta(mp.mpf(3)/4)*mp.dirichlet(mp.mpf(3)/4,[0,1,0,-1])/mp.zeta(mp.mpf(3)/2)
        Bsmall = -(4*mp.sqrt(2)/15)*Z/(2*mp.pi)
        return {'Zprim_regularized':float(Z), 'B_small_model':float(Bsmall),
                'B_large_mellin_model':float(Bsmall/4), 'B_total_model':float(5*Bsmall/4),
                'rho_integral_mellin_model':float(mp.pi*Z/5)}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    result = {'status':'MODEL_CALCULATION_NOT_HULL_PROOF', **model_coefficients(),
              'cutoffs':[{'T':T, 'integral':cutoff_integral(T)} for T in [50,100,250,500,1000,2000]],
              'rigorous_cutoff_error_bound':None}
    with args.out.open('x', encoding='utf-8') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
