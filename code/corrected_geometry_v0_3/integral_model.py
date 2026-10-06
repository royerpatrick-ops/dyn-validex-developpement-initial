"""Algorithmes de lignes du modele parabolique; aucune preuve du transport global.

G_cible estime uniquement la contribution de sa fenetre de cisaillement.
La couverture des autres configurations n'est pas certifiee.
"""
import math
import numpy as np

def validate_sampling(a, n, kmax):
    if not math.isfinite(a) or a <= 0:
        raise ValueError("a doit etre positif et fini")
    if not isinstance(n, (int, np.integer)) or isinstance(n, bool) or n < 2:
        raise ValueError("n doit etre un entier >= 2")
    if not isinstance(kmax, (int, np.integer)) or isinstance(kmax, bool) or kmax < 1:
        raise ValueError("kmax doit etre un entier >= 1")

def echantillon(a, n, rng, kmax=200000):
    """Pour un a donné : n tirages (δ, s, β) -> contribution c et indicatrice d'arête."""
    validate_sampling(a, n, kmax)
    delta = rng.uniform(0, 1 / a, n); s = rng.uniform(0, a, n); beta = rng.uniform(0, a, n)
    c = np.zeros(n); arete = np.zeros(n, bool)
    actif = np.ones(n, bool); k = 0
    while actif.any() and k < kmax:
        idx = np.where(actif)[0]
        d = delta[idx] + k / a
        r = np.sqrt(2 * d)
        off = np.mod(s[idx] + k * beta[idx], a)          # position du point le plus à gauche ≥ 0
        # points x = off + j a ; ceux dans [-r, r]
        jmin = np.ceil((-r - off) / a); jmax = np.floor((r - off) / a)
        npts = jmax - jmin + 1
        trouve = npts >= 1
        deux = npts >= 2
        sel = idx[trouve & deux]
        if sel.size:
            x1 = off[trouve & deux] + jmin[trouve & deux] * a
            x2 = off[trouve & deux] + jmax[trouve & deux] * a
            dd = d[trouve & deux]
            c[sel] = dd * (x2 - x1) - (x2 ** 3 - x1 ** 3) / 3
            arete[sel] = True
        actif[idx[trouve]] = False
        k += 1
    if actif.any():
        raise RuntimeError(f"a = {a}: {actif.sum()} tirages sans droite trouvée")
    return c, arete


def G_cible(a, n, rng, cb=5.0, cs=5.0, kmax=200000):
    validate_sampling(a, n, kmax)
    if not (math.isfinite(cb) and math.isfinite(cs) and cb > 0 and cs > 0):
        raise ValueError("fenetres positives et finies requises")
    if 2 * max(cb, cs) >= a ** 3:
        raise ValueError("fenetre recouvrant plusieurs periodes: poids invalide")
    Wb = cb / a ** 2; Ws = cs / a ** 2
    delta = rng.uniform(0, 1 / a, n)
    beta = rng.uniform(-Wb, Wb, n)
    K = (a * a / 8 - delta) * a                      # ligne où la corde atteint la largeur a
    s = np.mod(a / 2 - K * beta + rng.uniform(-Ws, Ws, n), a)
    poids = (2 * Wb / a) * (2 * Ws / a)
    c = np.zeros(n); act = np.arange(n); k = 0
    while act.size and k < kmax:
        d = delta[act] + k / a; r = np.sqrt(2 * d)
        off = np.mod(s[act] + k * beta[act], a)
        jmin = np.ceil((-r - off) / a); jmax = np.floor((r - off) / a); npt = jmax - jmin + 1
        tr = npt >= 1; dx = npt >= 2; sel = act[dx]
        x1 = off[dx] + jmin[dx] * a; x2 = off[dx] + jmax[dx] * a
        c[sel] = d[dx] * (x2 - x1) - (x2 ** 3 - x1 ** 3) / 3
        act = act[~tr]; k += 1
    if act.size:
        raise RuntimeError(f"a={a}: {act.size} tirages non resolus au kmax={kmax}")
    return c.mean() * poids, c.std(ddof=1) * poids / math.sqrt(n)
