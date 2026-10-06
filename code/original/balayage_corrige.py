"""
Balayage corrigé : déficit de périmètre de l'enveloppe convexe des points
entiers contenus dans une superellipse tournée |x/R|^p + |y/R|^p <= 1.

Corrections par rapport aux étapes 1E-1J :
  - aucune grille bornée : pour chaque colonne entière x, l'intervalle
    [y_min, y_max] est calculé exactement (bisection), donc la forme ne
    peut pas être tronquée ;
  - P_cont par intégration adaptative (scipy.quad) de la fonction d'appui ;
  - moyenne sur N_THETA orientations θ0 pour lisser les oscillations
    diophantiennes.
"""
import numpy as np  # pyright: ignore[reportMissingImports]
from scipy.spatial import ConvexHull  # pyright: ignore[reportMissingImports]
from scipy.integrate import quad  # pyright: ignore[reportMissingImports]
import time, sys

rng = np.random.default_rng(12345)
N_THETA = 16
THETAS = rng.uniform(0, np.pi / 2, N_THETA)

P_VALUES = [1.5, 2.0, 2.05, 2.1, 2.2, 2.5, 3.0, 4.0, 5.0, 10.0]
R_VALUES = np.unique(np.round(np.geomspace(30, 10000, 50)).astype(int)).astype(float)


def P_cont(R, p):
    """Formule de Cauchy : P = 2R ∫_0^π ||u(θ)||_q dθ, q = p/(p-1)."""
    q = p / (p - 1.0)
    f = lambda t: (abs(np.cos(t)) ** q + abs(np.sin(t)) ** q) ** (1.0 / q)
    v1, _ = quad(f, 0, np.pi / 2, limit=400, epsabs=1e-14, epsrel=1e-14)
    return 2 * R * 2 * v1  # symétrie [0,π/2] -> [0,π]


def F(x, y, R, p, c, s):
    return np.abs((c * x + s * y) / R) ** p + np.abs((-s * x + c * y) / R) ** p


def hull_points(R, p, th):
    c, s = np.cos(th), np.sin(th)
    ext = R * max(1.0, 2.0 ** (0.5 - 1.0 / p)) + 2
    xs = np.arange(-int(np.ceil(ext)), int(np.ceil(ext)) + 1, dtype=float)
    # 1) y* minimisant F(x,.) : recherche dorée (F convexe en y)
    lo = np.full_like(xs, -ext); hi = np.full_like(xs, ext)
    g = (np.sqrt(5) - 1) / 2
    for _ in range(80):
        a = hi - g * (hi - lo); b = lo + g * (hi - lo)
        fa = F(xs, a, R, p, c, s); fb = F(xs, b, R, p, c, s)
        m = fa < fb
        hi = np.where(m, b, hi); lo = np.where(m, lo, a)
    ystar = 0.5 * (lo + hi)
    inside = F(xs, ystar, R, p, c, s) <= 1.0
    xs, ystar = xs[inside], ystar[inside]
    # 2) bornes haute et basse par bisection
    def bound(direction):
        a = ystar.copy(); b = ystar + direction * 2 * ext
        for _ in range(70):
            mid = 0.5 * (a + b)
            ok = F(xs, mid, R, p, c, s) <= 1.0
            a = np.where(ok, mid, a); b = np.where(ok, b, mid)
        return a
    ytop = np.floor(bound(+1) + 1e-9)
    ybot = np.ceil(bound(-1) - 1e-9)
    # ajustement exact aux entiers (sécurité numérique)
    bad = F(xs, ytop, R, p, c, s) > 1.0; ytop[bad] -= 1
    bad = F(xs, ybot, R, p, c, s) > 1.0; ybot[bad] += 1
    ok = ytop >= ybot
    xs, ytop, ybot = xs[ok], ytop[ok], ybot[ok]
    return np.vstack([np.column_stack([xs, ytop]), np.column_stack([xs, ybot])])


def hull_perim(pts):
    h = ConvexHull(pts)
    v = pts[h.vertices]
    d = np.diff(np.vstack([v, v[:1]]), axis=0)
    return float(np.hypot(d[:, 0], d[:, 1]).sum()), len(h.vertices)


def main():
    rows = []
    t0 = time.time()
    for p in P_VALUES:
        for R in R_VALUES:
            Pc = P_cont(R, p)
            for th in THETAS:
                Ph, m = hull_perim(hull_points(R, p, th))
                assert Ph <= Pc + 1e-6, (p, R, th, Ph, Pc)
                rows.append((p, R, th, Pc, Ph, Pc - Ph, (Pc - Ph) / Pc, m))
        print(f"p={p:5.2f} fait ({time.time()-t0:.0f}s)", flush=True)
    arr = np.array(rows)
    np.savetxt("balayage_corrige.csv", arr, delimiter=",",
               header="p,R,theta0,P_cont,P_hull,deficit_abs,D_rel,m_sommets",
               comments="", fmt="%.10g")


if __name__ == "__main__":
    main()
