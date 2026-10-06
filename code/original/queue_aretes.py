"""Queue des arêtes longues : contribution à C_eff des arêtes de longueur normalisée ℓ' ≥ 2,
pour R de 10^3 à 10^6 (disque, réseau décalé aléatoire), et test de la coupure ℓ'_max ∝ R^(1/6)."""
import numpy as np, math
from decomposition_aretes import aretes
rng = np.random.default_rng(11)
plan = {1000: 4000, 10000: 1500, 100000: 400, 1000000: 120}
bins = 2.0 ** np.arange(1, 8, 0.5)
tab = {}
for R, N in plan.items():
    acc = np.zeros(len(bins) - 1); maxl = []
    for _ in range(N):
        th = rng.uniform(0, math.pi / 2); u, v = rng.uniform(0, 1, 2)
        L, c, gg = aretes(R, th, u, v)
        lp = L * R ** (-1 / 3)
        idx = np.digitize(lp, bins) - 1
        ok = (idx >= 0) & (idx < len(bins) - 1)
        np.add.at(acc, idx[ok], c[ok])
        maxl.append(lp.max())
    acc *= R ** (1 / 3) / (2 * np.pi) / N
    tab[R] = acc
    print(f"R = {R:>8d} ({N} tirages) : ℓ'max médian = {np.median(maxl):6.2f}"
          f"   ℓ'max / R^(1/6) = {np.median(maxl)/R**(1/6):.3f}")
np.save("queue_aretes.npy", {"bins": bins, **{str(k): v for k, v in tab.items()}}, allow_pickle=True)
print("\nPart de C_eff par demi-octave de ℓ' (lignes : bornes ; colonnes : R)")
print("   ℓ' ∈             " + "".join(f"{R:>11d}" for R in plan))
for k in range(len(bins) - 1):
    print(f"   [{bins[k]:6.2f},{bins[k+1]:6.2f})   " + "".join(f"{tab[R][k]:11.4f}" for R in plan))
