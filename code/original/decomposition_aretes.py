"""
Décomposition du déficit de périmètre arête par arête (disque de rayon R, réseau décalé).
Pour l'arête i entre deux sommets consécutifs, contribution c_i = R·Δφ_i − |e_i|,
où Δφ_i est l'angle (vu du centre) entre les deux sommets : Σ c_i = 2πR − P_hull = ΔP exactement.
On classe les arêtes par longueur normalisée ℓ' = |e| · κ_a^{1/3} = |e| R^{-1/3}
(échelle locale naturelle) et on somme les contributions par tranche dyadique de ℓ'.
Normalisation : C = ΔP · R^{1/3} / (2π), donc chaque tranche donne sa part de C_eff.
"""
import numpy as np, math, sys
from scipy.spatial import ConvexHull
from balayage_translations import columns

def aretes(R, th, u, v):
    pts = columns("ellipse", 1.0, float(R), th, u, v)
    h = ConvexHull(pts); w = pts[h.vertices]            # sommets, sens trigonométrique
    e = np.roll(w, -1, axis=0) - w
    L = np.hypot(e[:, 0], e[:, 1])
    phi = np.unwrap(np.arctan2(w[:, 1], w[:, 0]))
    dphi = np.diff(np.append(phi, phi[0] + 2 * np.pi))
    c = R * dphi - L
    # vecteur primitif et multiplicité de chaque arête (coordonnées entières : sommets = k - (u,v))
    ei = np.rint(e).astype(np.int64)
    g = np.gcd(np.abs(ei[:, 0]), np.abs(ei[:, 1]))
    return L, c, g

def main():
    bins = 2.0 ** np.arange(-6, 12)
    rng = np.random.default_rng(5)
    print("Part de C_eff portée par chaque tranche de longueur d'arête ℓ' = |e|·R^(-1/3)")
    res = {}
    for R in [10**4, 10**5, 10**6]:
        N = 64 if R < 10**6 else 32
        acc = np.zeros(len(bins) - 1); nb = np.zeros(len(bins) - 1); tot = 0; mult = []
        for _ in range(N):
            th = rng.uniform(0, math.pi / 2); u, v = rng.uniform(0, 1, 2)
            L, c, g = aretes(R, th, u, v)
            lp = L * R ** (-1 / 3)
            idx = np.digitize(lp, bins) - 1
            for k in range(len(bins) - 1):
                m = idx == k; acc[k] += c[m].sum(); nb[k] += m.sum()
            tot += c.sum(); mult.append(np.mean(g > 1))
        norm = R ** (1 / 3) / (2 * np.pi) / N
        res[R] = (acc * norm, nb / N, tot * norm)
        print(f"\nR = {R:>8d}   C_eff = {tot*norm:.4f}   arêtes non primitives : {100*np.mean(mult):.1f} %")
        print("   ℓ' ∈          part de C     arêtes/échantillon")
        for k in range(len(bins) - 1):
            if nb[k] > 0:
                print(f"   [{bins[k]:7.3f},{bins[k+1]:7.3f})  {acc[k]*norm:8.4f}   {nb[k]/N:10.1f}")
    np.save("decomposition_aretes.npy", {R: res[R] for R in res}, allow_pickle=True)


if __name__ == "__main__":
    main()
