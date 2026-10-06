"""Contribution à C_eff des arêtes du disque réel, classées par leur vecteur primitif w0 ∈ Z²
(longueur r0), comparée à la prédiction du régime « petits a » : n(r0)·K·R^(-1/6)·r0^(-3/2)/(2π)."""
import math, numpy as np
from collections import defaultdict
from decomposition_aretes import aretes
from scipy.spatial import ConvexHull
from balayage_translations import columns
K = 0.629515
def classes(R, N, rng):
    acc = defaultdict(float)
    for _ in range(N):
        th = rng.uniform(0, math.pi/2); u, v = rng.uniform(0, 1, 2)
        pts = columns("ellipse", 1.0, float(R), th, u, v)
        h = ConvexHull(pts); w = pts[h.vertices]
        e = np.roll(w, -1, axis=0) - w; L = np.hypot(e[:,0], e[:,1])
        phi = np.unwrap(np.arctan2(w[:,1], w[:,0])); dphi = np.diff(np.append(phi, phi[0]+2*np.pi))
        c = R*dphi - L
        ei = np.rint(e).astype(np.int64); g = np.gcd(np.abs(ei[:,0]), np.abs(ei[:,1]))
        r2 = (ei[:,0]//g)**2 + (ei[:,1]//g)**2
        for k, ck in zip(r2, c): acc[int(k)] += ck
    return {k: v*R**(1/3)/(2*math.pi)/N for k, v in acc.items()}
nprim = defaultdict(int)
for p in range(-30, 31):
    for q in range(-30, 31):
        if (p or q) and math.gcd(p, q) == 1: nprim[p*p+q*q] += 1
rng = np.random.default_rng(4)
for R, N in [(10**5, 200), (10**6, 60)]:
    m = classes(R, N, rng)
    print(f"\nR = {R}  (ε^1/2 = R^-1/6 = {R**(-1/6):.4f})")
    print("  r0²   n   mesuré     prédit    rapport")
    for k in sorted(nprim)[:10]:
        pred = nprim[k]*K*R**(-1/6)*k**(-0.75)/(2*math.pi)
        print(f"  {k:3d} {nprim[k]:3d}  {m.get(k,0):.5f}   {pred:.5f}   {m.get(k,0)/pred:.3f}")
