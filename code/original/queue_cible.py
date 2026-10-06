"""Contribution des grands a : échantillonnage dans la fenêtre β ≈ 0, s + Kβ ≈ a/2."""
import math, sys
import numpy as np  # pyright: ignore[reportMissingImports]

def G_cible(a, n, rng, cb=10.0, cs=10.0):
    Wb = cb / a ** 2; Ws = cs / a ** 2
    delta = rng.uniform(0, 1 / a, n)
    beta = rng.uniform(-Wb, Wb, n)
    K = (a * a / 8 - delta) * a                      # ligne où la corde atteint la largeur a
    s = np.mod(a / 2 - K * beta + rng.uniform(-Ws, Ws, n), a)
    poids = (2 * Wb / a) * (2 * Ws / a)
    c = np.zeros(n); act = np.arange(n); k = 0
    while act.size:
        d = delta[act] + k / a; r = np.sqrt(2 * d)
        off = np.mod(s[act] + k * beta[act], a)
        jmin = np.ceil((-r - off) / a); jmax = np.floor((r - off) / a); npt = jmax - jmin + 1
        tr = npt >= 1; dx = npt >= 2; sel = act[dx]
        x1 = off[dx] + jmin[dx] * a; x2 = off[dx] + jmax[dx] * a
        c[sel] = d[dx] * (x2 - x1) - (x2 ** 3 - x1 ** 3) / 3
        act = act[~tr]; k += 1
    return c.mean() * poids, c.std() * poids / math.sqrt(n)

def main():
    rng = np.random.default_rng(2)
    print("contrôle de fenêtre (a = 10) :")
    for cb, cs in [(5, 5), (10, 10), (20, 20), (40, 40)]:
        g, e = G_cible(10.0, 400000, rng, cb, cs); print(f"   fenêtres ×{cb:>2}: a³·G = {1000*g:.4f} ± {1000*e:.4f}")
    print("a³·G(a) en fonction de a :")
    for a in [6, 8, 10, 14, 20, 28, 40, 56, 80]:
        g, e = G_cible(float(a), 300000, rng, 20, 20); print(f"   a = {a:3d} : a³·G = {a**3*g:.4f} ± {a**3*e:.4f}")


if __name__ == "__main__":
    main()
