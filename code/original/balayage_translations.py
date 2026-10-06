"""
Balayage v2 : déficit de périmètre de l'enveloppe entière avec moyenne
sur rotations ET translations aléatoires du centre.

Correction par rapport à balayage_corrige.py : pour le disque (et les formes
presque circulaires), une rotation autour d'un centre entier ne change pas
l'ensemble des points du réseau ; la moyenne sur les rotations seules ne
moyenne donc rien. On tire ici, pour chaque échantillon, une rotation
θ ∈ [0, π/2) ET un décalage (u, v) ∈ [0, 1)² du réseau.

Sorties : balayage_translations.csv
  forme, param, R, echantillon, theta, u, v, P_cont, P_hull, deficit_abs, m_sommets
"""
import sys, time, math, zlib
import numpy as np  # pyright: ignore[reportMissingImports]
from scipy.spatial import ConvexHull  # pyright: ignore[reportMissingImports]
from scipy.integrate import quad  # pyright: ignore[reportMissingImports]
from multiprocessing import Pool

N_SAMPLES = 64
R_VALUES = np.unique(np.round(np.geomspace(100, 30000, 26)).astype(int)).astype(float)
SHAPES = ([("super", p) for p in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]]
          + [("ellipse", k) for k in [1.5, 2.0, 3.0, 5.0]])


# ------------------------------------------------------------------
# Invariants exacts I_a = ∫ κ^a ds  (superellipse : paramétrage par x,
# changement x = s^6 qui supprime les singularités de courbure)
# ------------------------------------------------------------------
def I_super(p, a):
    def f(x):
        if x <= 0.0:
            return 0.0 if a > 0 else 1.0
        w = 1.0 - x ** p
        yp = -x ** (p - 1) * w ** (1 / p - 1)
        ypp = -(p - 1) * x ** (p - 2) * w ** (1 / p - 2)
        g = 1.0 + yp * yp
        return (abs(ypp) / g ** 1.5) ** a * math.sqrt(g)
    xm = 2.0 ** (-1.0 / p)
    v, _ = quad(lambda s: f(s ** 6) * 6 * s ** 5, 0.0, xm ** (1 / 6),
                limit=2000, epsabs=1e-14, epsrel=1e-13)
    return 8.0 * v


def I_ellipse(k, a):
    A, B = 1.0, 1.0 / k
    sp = lambda t: math.hypot(A * math.sin(t), B * math.cos(t))
    v, _ = quad(lambda t: (A * B / sp(t) ** 3) ** a * sp(t), 0, math.pi / 2,
                limit=1000, epsabs=1e-14, epsrel=1e-13)
    return 4.0 * v


def invariants(kind, prm):
    I = I_super if kind == "super" else I_ellipse
    return I(prm, 0.0), I(prm, 1 / 3), I(prm, 4 / 3)  # P, Omega, J


# ------------------------------------------------------------------
# Points extrêmes de chaque colonne du réseau Z² - (u, v) dans R·K
# ------------------------------------------------------------------
def columns(kind, prm, R, th, u, v):
    c, s = math.cos(th), math.sin(th)
    if kind == "super":
        p = prm
        F = lambda x, y: (np.abs((c * x + s * y) / R) ** p
                          + np.abs((-s * x + c * y) / R) ** p)
        ext = R * max(1.0, 2.0 ** (0.5 - 1.0 / p)) + 2
    else:
        a, b = R, R / prm
        A = c * c / a ** 2 + s * s / b ** 2
        Bc = c * s * (1 / a ** 2 - 1 / b ** 2)
        C = s * s / a ** 2 + c * c / b ** 2
        F = lambda x, y: A * x * x + 2 * Bc * x * y + C * y * y
        ext = math.sqrt(C / (A * C - Bc * Bc)) + 2
    kx = np.arange(math.floor(-ext + u), math.ceil(ext + u) + 1, dtype=float)
    xs = kx - u
    if kind == "super":
        lo = np.full_like(xs, -ext); hi = np.full_like(xs, ext)
        g = (math.sqrt(5) - 1) / 2
        for _ in range(90):
            a1 = hi - g * (hi - lo); b1 = lo + g * (hi - lo)
            m = F(xs, a1) < F(xs, b1)
            hi = np.where(m, b1, hi); lo = np.where(m, lo, a1)
        ym = 0.5 * (lo + hi)
        ok = F(xs, ym) <= 1.0
        xs, ym = xs[ok], ym[ok]

        def bound(sign):
            a1 = ym.copy(); b1 = ym + sign * 2 * ext
            for _ in range(75):
                mid = 0.5 * (a1 + b1)
                inside = F(xs, mid) <= 1.0
                a1 = np.where(inside, mid, a1); b1 = np.where(inside, b1, mid)
            return a1
        ytop, ybot = bound(+1), bound(-1)
    else:
        disc = Bc * Bc * xs * xs - C * (A * xs * xs - 1)
        ok = disc >= 0
        xs, disc = xs[ok], disc[ok]
        ytop = (-Bc * xs + np.sqrt(disc)) / C
        ybot = (-Bc * xs - np.sqrt(disc)) / C
    # entiers k_y tels que k_y - v soit dans [ybot, ytop], test exact
    kt = np.floor(ytop + v + 1e-9); kb = np.ceil(ybot + v - 1e-9)
    for _ in range(2):
        kt = np.where(F(xs, kt - v) > 1.0, kt - 1, kt)
        kb = np.where(F(xs, kb - v) > 1.0, kb + 1, kb)
    keep = kt >= kb
    xs, kt, kb = xs[keep], kt[keep], kb[keep]
    return np.vstack([np.column_stack([xs, kt - v]), np.column_stack([xs, kb - v])])


def hull_perim(pts):
    h = ConvexHull(pts)
    w = pts[h.vertices]
    d = np.diff(np.vstack([w, w[:1]]), axis=0)
    return float(np.hypot(d[:, 0], d[:, 1]).sum()), len(h.vertices)


# ------------------------------------------------------------------
def job(args):
    kind, prm, R, P1 = args
    rng = np.random.default_rng(zlib.crc32(f"{kind}-{prm}-{int(R)}".encode()))
    out = []
    for i in range(N_SAMPLES):
        th = rng.uniform(0, math.pi / 2); u, v = rng.uniform(0, 1, 2)
        Ph, m = hull_perim(columns(kind, prm, R, th, u, v))
        Pc = P1 * R
        assert Ph <= Pc + 1e-6
        out.append((kind, prm, R, i, th, u, v, Pc, Ph, Pc - Ph, m))
    return out


def main():
    inv = {sh: invariants(*sh) for sh in SHAPES}
    tasks = [(k, p, R, inv[(k, p)][0]) for (k, p) in SHAPES for R in R_VALUES]
    tasks.sort(key=lambda t: -t[2])  # gros R d'abord (équilibrage)
    t0 = time.time(); rows = []
    with Pool() as pool:
        for n, res in enumerate(pool.imap_unordered(job, tasks), 1):
            rows.extend(res)
            if n % 20 == 0:
                print(f"{n}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)
    with open("balayage_translations.csv", "w") as f:
        f.write("forme,param,R,echantillon,theta,u,v,P_cont,P_hull,deficit_abs,m_sommets\n")
        for r in sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3])):
            f.write(",".join(str(x) for x in r) + "\n")
    with open("invariants.csv", "w") as f:
        f.write("forme,param,P,Omega,I43\n")
        for (k, p), (P, Om, J) in inv.items():
            f.write(f"{k},{p},{P:.12f},{Om:.12f},{J:.12f}\n")
    print(f"Terminé en {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
