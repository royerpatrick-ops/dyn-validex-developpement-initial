"""
enveloppe_commune.py — routines communes aux scripts de validation (version corrigée,
3 octobre 2026).

Ce module remplace la fonction de bornes de colonnes des anciens scripts.

Bug corrigé
-----------
L'ancienne fonction cherchait un point de départ dans chaque colonne en testant
600 valeurs de y réparties sur [-R√2-2, R√2+2], soit un pas d'environ 0,0047·R
(472 mailles à R = 10⁵). Aux extrémités gauche et droite de la forme, une colonne
dont la corde est plus courte que ce pas ne contient aucun point de la grille :
elle était déclarée vide et disparaissait de l'enveloppe. Le périmètre de
l'enveloppe était alors sous-estimé, d'où un biais positif sur C_eff, croissant
avec R, fort pour une ellipse tournée et nul pour le disque (dont les colonnes
extrêmes passent par y = 0, valeur ajoutée explicitement à la grille).

Seconde correction : le test d'appartenance utilisait F ≤ 1 + 1e-9, ce qui
agrandit la forme d'environ 5·10⁻¹⁰·R et biaise C_eff vers le bas (−0,002 sur le
disque à R = 10⁵). Le test est désormais strict : F ≤ niveau.

Méthode
-------
* Forme convexe : pour chaque colonne, minimum de F le long de la colonne par
  section dorée (F est unimodale le long d'une droite pour toutes les formes
  convexes utilisées), test direct des entiers voisins de ce minimum, puis
  bisection ENTIÈRE vers le haut et vers le bas. Le résultat est exact par
  rapport au test d'appartenance : aucune tolérance, aucune grille.
* Forme non convexe (rosace, lemniscate) : énumération de tous les points du
  réseau, colonne par colonne. Coût en R², réservé aux R ≤ 10⁴.

Convention de position (celle des scripts d'origine) : le point (i, j) du
réseau est dans RK si F((i cos θ + j sin θ − u)/R, (−i sin θ + j cos θ − v)/R)
≤ niveau.

Dépendances : numpy, scipy.
"""
import math
import numpy as np
from scipy.integrate import quad
from scipy.spatial import ConvexHull

C_INF = 0.7192
B_CONST = 0.2894


def C_theo(R, ratio):
    """Loi à deux termes : C_eff(R) = C∞ − B · (I_3/2 / I_4/3) · R^(−1/6)."""
    return C_INF - B_CONST * float(ratio) * float(R) ** (-1.0 / 6.0)


# ============================================================================
# 1. FORMES
# ============================================================================

class Forme:
    """Une forme du plan, à l'échelle R = 1.

    F(X, Y)        : fonction dont le sous-niveau {F ≤ niveau} est la forme
    derivees(X, Y) : (Fx, Fy, Fxx, Fxy, Fyy), pour la courbure du bord
    courbe(n)      : n points du bord, parcourus dans l'ordre
    rho_max        : rayon d'un disque centré à l'origine contenant la forme
    convexe        : True si la forme est convexe
    """

    def __init__(self, nom, F, niveau, derivees, courbe, rho_max, convexe):
        self.nom, self.F, self.niveau = nom, F, niveau
        self.derivees, self.courbe = derivees, courbe
        self.rho_max, self.convexe = rho_max, convexe

    def dedans(self, X, Y):
        return self.F(X, Y) <= self.niveau

    def courbure(self, X, Y):
        Fx, Fy, Fxx, Fxy, Fyy = self.derivees(X, Y)
        return np.abs(Fxx * Fy**2 - 2 * Fxy * Fx * Fy + Fyy * Fx**2) / (Fx**2 + Fy**2) ** 1.5


def _lame_generique(nom, a, b, p):
    """|X/a|^p + |Y/b|^p ≤ 1 (superellipse, ellipse, courbe de Lamé)."""
    ap, bp = a**p, b**p

    def F(X, Y):
        return np.abs(X) ** p / ap + np.abs(Y) ** p / bp

    def derivees(X, Y):
        ax, ay = np.abs(X), np.abs(Y)
        Fx = p * np.sign(X) * ax ** (p - 1) / ap
        Fy = p * np.sign(Y) * ay ** (p - 1) / bp
        with np.errstate(divide="ignore", invalid="ignore"):
            Fxx = p * (p - 1) * ax ** (p - 2) / ap
            Fyy = p * (p - 1) * ay ** (p - 2) / bp
        return Fx, Fy, Fxx, np.zeros_like(X), Fyy

    def courbe(n):
        t = (np.arange(n) + 0.5) * 2 * np.pi / n
        c, s = np.cos(t), np.sin(t)
        return a * np.sign(c) * np.abs(c) ** (2 / p), b * np.sign(s) * np.abs(s) ** (2 / p)

    rho = math.hypot(a, b) if p >= 2 else max(a, b)   # borne simple et sûre
    return Forme(nom, F, 1.0, derivees, courbe, rho, convexe=(p >= 1))


def superellipse(p):
    f = _lame_generique(f"superellipse p={p}", 1.0, 1.0, p)
    f.p = p
    return f


def ellipse(ratio):
    """Ellipse de demi-axes (1, 1/ratio), comme dans les scripts d'origine."""
    f = _lame_generique(f"ellipse a/b={ratio}", 1.0, 1.0 / ratio, 2.0)
    f.ratio = ratio
    return f


def lame(a, b, p):
    return _lame_generique(f"Lamé a={a} b={b} p={p}", a, b, p)


def cassini(a, c):
    """Ovale de Cassini (X²+Y²)² − 2c²(X²−Y²) ≤ a⁴ − c⁴ ; convexe si c/a ≤ 1/√2."""
    convexe = c / a <= 1 / math.sqrt(2) + 1e-12

    def F(X, Y):
        r2 = X * X + Y * Y
        return r2 * r2 - 2 * c * c * (X * X - Y * Y)

    def derivees(X, Y):
        r2 = X * X + Y * Y
        return (4 * X * r2 - 4 * c * c * X, 4 * Y * r2 + 4 * c * c * Y,
                4 * r2 + 8 * X * X - 4 * c * c, 8 * X * Y, 4 * r2 + 8 * Y * Y + 4 * c * c)

    def courbe(n):
        th = (np.arange(n) + 0.5) * 2 * np.pi / n
        u = c**2 * np.cos(2 * th) + np.sqrt(c**4 * np.cos(2 * th) ** 2 + a**4 - c**4)
        r = np.sqrt(u)
        return r * np.cos(th), r * np.sin(th)

    return Forme(f"Cassini a={a} c={c}", F, a**4 - c**4, derivees, courbe,
                 math.sqrt(a * a + c * c), convexe)


def lemniscate():
    """Lemniscate de Bernoulli (X²+Y²)² ≤ X² − Y² (non convexe, deux lobes)."""

    def F(X, Y):
        r2 = X * X + Y * Y
        return r2 * r2 - (X * X - Y * Y)

    def derivees(X, Y):
        r2 = X * X + Y * Y
        return (4 * X * r2 - 2 * X, 4 * Y * r2 + 2 * Y,
                4 * r2 + 8 * X * X - 2, 8 * X * Y, 4 * r2 + 8 * Y * Y + 2)

    def courbe(n):
        m = n // 2
        th = -np.pi / 4 + (np.arange(m) + 0.5) * (np.pi / 2) / m
        th = np.concatenate([th, th + np.pi])
        r = np.sqrt(np.maximum(np.cos(2 * th), 0))
        return r * np.cos(th), r * np.sin(th)

    return Forme("lemniscate", F, 0.0, derivees, courbe, 1.0, convexe=False)


def rosace(k):
    """Rosace r ≤ |cos kθ|, soit (X²+Y²)^((k+1)/2) ≤ |Re (X+iY)^k| (non convexe)."""
    m = (k + 1) / 2

    def F(X, Y):
        z = X + 1j * Y
        return (X * X + Y * Y) ** m - np.abs((z**k).real)

    def derivees(X, Y):
        z = X + 1j * Y
        r2 = X * X + Y * Y
        s = np.sign((z**k).real)
        d1 = k * z ** (k - 1)
        d2 = k * (k - 1) * z ** (k - 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            Fx = 2 * m * X * r2 ** (m - 1) - s * d1.real
            Fy = 2 * m * Y * r2 ** (m - 1) + s * d1.imag
            Fxx = 2 * m * r2 ** (m - 1) + 4 * m * (m - 1) * X * X * r2 ** (m - 2) - s * d2.real
            Fyy = 2 * m * r2 ** (m - 1) + 4 * m * (m - 1) * Y * Y * r2 ** (m - 2) + s * d2.real
            Fxy = 4 * m * (m - 1) * X * Y * r2 ** (m - 2) + s * d2.imag
        return Fx, Fy, Fxx, Fxy, Fyy

    def courbe(n):
        th = (np.arange(n) + 0.5) * 2 * np.pi / n
        r = np.abs(np.cos(k * th))
        return r * np.cos(th), r * np.sin(th)

    return Forme(f"rosace k={k}", F, 0.0, derivees, courbe, 1.0, convexe=False)


# ============================================================================
# 2. INVARIANTS ET PÉRIMÈTRE CONTINU (à R = 1)
# ============================================================================

def invariants_superellipse(p):
    """(P, I_4/3, I_3/2) de |x|^p + |y|^p ≤ 1, par quadrature sur un huitième du bord.

    Changement de variable x = w⁶ pour absorber les singularités de courbure
    (courbure infinie pour p < 2, nulle pour p > 2, aux points des axes).
    """
    if abs(p - 2) < 1e-12:
        return 2 * math.pi, 2 * math.pi, 2 * math.pi
    t = 2 ** (-1 / p)

    def I(e):
        def g(w):
            x = w**6
            if x <= 0:
                return 0.0
            fp = x ** (p - 1) * (1 - x**p) ** (1 / p - 1)
            fpp = (p - 1) * x ** (p - 2) * (1 - x**p) ** (1 / p - 2)
            return fpp**e * (1 + fp * fp) ** ((1 - 3 * e) / 2) * 6 * w**5
        return 8 * quad(g, 0, t ** (1 / 6), epsabs=0, epsrel=1e-12, limit=400)[0]

    return I(0.0), I(4 / 3), I(1.5)


def invariants_ellipse(ratio):
    """(P, I_4/3, I_3/2) de l'ellipse de demi-axes (1, 1/ratio)."""
    a, b = 1.0, 1.0 / ratio

    def I(e):
        f = lambda t: (a * b) ** e * (a * a * math.sin(t) ** 2 + b * b * math.cos(t) ** 2) ** ((1 - 3 * e) / 2)
        return 4 * quad(f, 0, math.pi / 2, epsabs=0, epsrel=1e-13, limit=200)[0]

    return I(0.0), I(4 / 3), I(1.5)


def invariants_bord(forme, n=2_000_000):
    """(P_conv, I_4/3, I_3/2) de l'enveloppe convexe continue de la forme, à R = 1.

    Le bord est échantillonné en n points ; la courbure est calculée exactement en
    chaque point par la formule implicite (dérivées de F), et ds par les cordes.
    Pour une forme non convexe, seuls les arcs de courbe qui appartiennent au bord
    de l'enveloppe convexe comptent dans les intégrales ; les segments droits de
    l'enveloppe ont une courbure nulle. Convient aux bords sans point de courbure
    infinie (pour les superellipses, utiliser invariants_superellipse).
    """
    x, y = forme.courbe(n)
    nxt = np.roll(np.arange(n), -1)
    seg = np.hypot(x[nxt] - x, y[nxt] - y)
    if forme.convexe:
        sur_bord = np.ones(n, dtype=bool)
        P = math.fsum(seg.tolist())
    else:
        pts = np.column_stack([x, y])
        hull = ConvexHull(pts)
        sur_bord = np.zeros(n, dtype=bool)
        sur_bord[hull.vertices] = True
        v = pts[hull.vertices]
        d = np.roll(v, -1, axis=0) - v
        P = math.fsum(np.hypot(d[:, 0], d[:, 1]).tolist())
    arc = sur_bord & sur_bord[nxt]
    k = np.zeros(n)
    k[sur_bord] = forme.courbure(x[sur_bord], y[sur_bord])

    def I(e):
        ke = k**e
        return math.fsum((0.5 * (ke + ke[nxt]) * seg)[arc].tolist())

    return P, I(4 / 3), I(1.5)


# ============================================================================
# 3. BORNES ENTIÈRES PAR COLONNE
# ============================================================================

_OR = (math.sqrt(5) - 1) / 2


def bornes_colonnes(xs, R, forme, theta, u, v, n_or=72):
    """Bornes entières (y_bas, y_haut) des points du réseau dans RK, colonne par colonne.

    xs : abscisses entières des colonnes (tableau de float, comme dans les scripts).
    Retourne deux tableaux de float, NaN pour les colonnes vides.
    """
    if not forme.convexe:
        return bornes_colonnes_brutes(xs, R, forme, theta, u, v)

    c, s = math.cos(theta), math.sin(theta)
    R = float(R)
    xs = np.asarray(xs, dtype=np.float64)
    Ymax = forme.rho_max * R + 3.0
    ntot = len(xs)

    def Fcol(x, y):
        return forme.F((c * x + s * y - u) / R, (-s * x + c * y - v) / R)

    # 3.1 Minimum de F le long de chaque colonne (section dorée vectorisée)
    a = np.full(ntot, -Ymax); b = np.full(ntot, Ymax)
    yc = b - _OR * (b - a); yd = a + _OR * (b - a)
    fc = Fcol(xs, yc); fd = Fcol(xs, yd)
    for _ in range(n_or):
        m = fc < fd
        a = np.where(m, a, yc)
        b = np.where(m, yd, b)
        y_new = np.where(m, b - _OR * (b - a), a + _OR * (b - a))
        f_new = Fcol(xs, y_new)
        yc, yd, fc, fd = (np.where(m, y_new, yd), np.where(m, yc, y_new),
                          np.where(m, f_new, fd), np.where(m, fc, f_new))
    y_etoile = 0.5 * (a + b)

    # 3.2 Graine entière : un des entiers voisins du minimum est dans RK
    #     si et seulement si la colonne contient un point du réseau.
    base = np.floor(y_etoile)
    graine = np.full(ntot, np.nan)
    for dj in (0.0, 1.0, -1.0, 2.0):
        cand = base + dj
        ok = np.isnan(graine) & forme.dedans((c * xs + s * cand - u) / R, (-s * xs + c * cand - v) / R)
        graine = np.where(ok, cand, graine)
    valide = ~np.isnan(graine)

    y_bas = np.full(ntot, np.nan); y_haut = np.full(ntot, np.nan)
    if not valide.any():
        return y_bas, y_haut
    xv = xs[valide]
    j0 = graine[valide].astype(np.int64)
    portee = int(math.ceil(2 * Ymax)) + 4

    def dedans_col(jj):
        jj = jj.astype(np.float64)
        return forme.dedans((c * xv + s * jj - u) / R, (-s * xv + c * jj - v) / R)

    # 3.3 Bisection entière vers le haut : lo dans RK, hi hors de RK
    lo, hi = j0.copy(), j0 + portee
    while np.any(hi - lo > 1):
        mid = (lo + hi) // 2
        ins = dedans_col(mid)
        lo = np.where(ins, mid, lo); hi = np.where(ins, hi, mid)
    y_haut[valide] = lo
    # 3.4 Bisection entière vers le bas : hi dans RK, lo hors de RK
    lo, hi = j0 - portee, j0.copy()
    while np.any(hi - lo > 1):
        mid = (lo + hi) // 2
        ins = dedans_col(mid)
        hi = np.where(ins, mid, hi); lo = np.where(ins, lo, mid)
    y_bas[valide] = hi
    return y_bas, y_haut


def bornes_colonnes_brutes(xs, R, forme, theta, u, v, max_elements=4_000_000):
    """Énumération complète (formes non convexes, ou contrôle) : coût en R²."""
    c, s = math.cos(theta), math.sin(theta)
    R = float(R)
    xs = np.asarray(xs, dtype=np.float64)
    Ymax = int(math.ceil(forme.rho_max * R + 3.0))
    ys = np.arange(-Ymax, Ymax + 1, dtype=np.float64)
    y_bas = np.full(len(xs), np.nan); y_haut = np.full(len(xs), np.nan)
    pas = max(1, max_elements // len(ys))
    for i0 in range(0, len(xs), pas):
        X = xs[i0:i0 + pas, None]
        M = forme.dedans((c * X + s * ys[None, :] - u) / R, (-s * X + c * ys[None, :] - v) / R)
        a_un = M.any(axis=1)
        prem = np.argmax(M, axis=1)
        dern = len(ys) - 1 - np.argmax(M[:, ::-1], axis=1)
        y_bas[i0:i0 + pas] = np.where(a_un, ys[prem], np.nan)
        y_haut[i0:i0 + pas] = np.where(a_un, ys[dern], np.nan)
    return y_bas, y_haut


def colonnes(R, forme):
    """Abscisses des colonnes à examiner (couvre toute position de la forme)."""
    xb = int(math.ceil(forme.rho_max * R)) + 3
    return np.arange(-xb, xb + 1, dtype=np.float64)


# ============================================================================
# 4. PÉRIMÈTRE DE L'ENVELOPPE
# ============================================================================

def hull_perimeter_from_bounds(xs, y_bot, y_top):
    """Périmètre de l'enveloppe convexe des extrémités de colonnes (inchangé, avec fsum)."""
    ok = ~np.isnan(y_top) & ~np.isnan(y_bot)
    xv = xs[ok]
    if len(xv) == 0:
        return 0.0
    pts = np.unique(np.vstack([np.column_stack([xv, y_top[ok]]),
                               np.column_stack([xv, y_bot[ok]])]), axis=0)
    if len(pts) < 3:
        if len(pts) == 2:
            return 2.0 * float(np.hypot(*(pts[1] - pts[0])))
        return 0.0
    try:
        v = pts[ConvexHull(pts).vertices]
    except Exception:            # points alignés
        d = pts.max(axis=0) - pts.min(axis=0)
        return 2.0 * float(np.hypot(*d))
    d = np.roll(v, -1, axis=0) - v
    return math.fsum(np.hypot(d[:, 0], d[:, 1]).tolist())


def perimetre_tirage(R, forme, theta, u, v, xs=None, n_or=72):
    """Périmètre de l'enveloppe entière pour une position (θ, u, v).

    n_or : itérations de section dorée. 48 suffisent (précision ≈ 10⁻⁴ maille à
    R = 10⁶, alors qu'une précision de 1 suffit pour la graine) ; les bornes
    finales sont exactes quel que soit n_or, par la bisection entière.
    """
    if xs is None:
        xs = colonnes(R, forme)
    yb, yt = bornes_colonnes(xs, R, forme, theta, u, v, n_or=n_or)
    return hull_perimeter_from_bounds(xs, yb, yt)


# ============================================================================
# 5. CONTRÔLES INTÉGRÉS
# ============================================================================

def _bornes_ellipse_fermees(xs, R, ratio, theta, u, v):
    """Bornes exactes de l'ellipse par les racines du trinôme en j (contrôle)."""
    c, s = math.cos(theta), math.sin(theta)
    A = c * xs - u; Bv = -s * xs - v
    al = s * s + ratio**2 * c * c
    be = s * A + ratio**2 * c * Bv
    ga = A * A + ratio**2 * Bv * Bv - R * R
    D = be * be - al * ga
    ok = D >= 0
    sq = np.sqrt(np.where(ok, D, 0))
    f = ellipse(ratio)
    ins = lambda jj: f.dedans((c * xs + s * jj - u) / R, (-s * xs + c * jj - v) / R)
    jlo = np.ceil((-be - sq) / al); jhi = np.floor((-be + sq) / al)
    for _ in range(2):
        jlo = np.where(ins(jlo - 1), jlo - 1, jlo); jlo = np.where(~ins(jlo), jlo + 1, jlo)
        jhi = np.where(ins(jhi + 1), jhi + 1, jhi); jhi = np.where(~ins(jhi), jhi - 1, jhi)
    ok &= jlo <= jhi
    return np.where(ok, jlo, np.nan), np.where(ok, jhi, np.nan)


def controle(verbeux=True):
    """Deux contrôles des bornes de colonnes. Lève une erreur en cas d'échec.

    1. Petits R : identité colonne par colonne avec l'énumération complète,
       pour des formes convexes variées.
    2. Grand R (10⁵) : identité colonne par colonne avec les bornes en forme
       fermée de l'ellipse 2:1, y compris pour un grand axe presque horizontal,
       le cas où l'ancienne grille perdait des colonnes.
    """
    rng = np.random.default_rng(2026)
    formes = [superellipse(1.5), superellipse(3.0), superellipse(10.0), ellipse(2.0),
              ellipse(5.0), lame(1.0, 0.5, 3.0), cassini(1.0, 0.7)]
    n1 = 0
    for k in range(42):
        f = formes[k % len(formes)]
        R = float(rng.uniform(20, 250)); th = rng.uniform(0, np.pi / 2); u, v = rng.uniform(0, 1, 2)
        xs = colonnes(R, f)
        a = bornes_colonnes(xs, R, f, th, u, v)
        b = bornes_colonnes_brutes(xs, R, f, th, u, v)
        if not (np.array_equal(a[0], b[0], equal_nan=True) and np.array_equal(a[1], b[1], equal_nan=True)):
            raise RuntimeError(f"Contrôle 1 en échec : {f.nom}, R={R:.1f}, θ={th:.3f}")
        n1 += 1
    n2 = 0
    for th in (0.002, 0.05, 0.4, 1.3):
        u, v = rng.uniform(0, 1, 2)
        R = 1e5; f = ellipse(2.0); xs = colonnes(R, f)
        a = bornes_colonnes(xs, R, f, th, u, v)
        b = _bornes_ellipse_fermees(xs, R, 2.0, th, u, v)
        if not (np.array_equal(a[0], b[0], equal_nan=True) and np.array_equal(a[1], b[1], equal_nan=True)):
            raise RuntimeError(f"Contrôle 2 en échec : ellipse 2:1, R=1e5, θ={th}")
        n2 += 1
    if verbeux:
        print(f"Contrôles des bornes : {n1} cas à petit R identiques à l'énumération complète, "
              f"{n2} cas à R = 10⁵ identiques aux bornes exactes de l'ellipse. OK\n")
    return True


if __name__ == "__main__":
    controle()
    for p in (1.5, 3.0, 10.0):
        P, I43, I32 = invariants_superellipse(p)
        print(f"superellipse p={p:<4}  P={P:.6f}  I43={I43:.4f}  I32={I32:.4f}")
    for r in (1.5, 2.0, 5.0):
        P, I43, I32 = invariants_ellipse(r)
        print(f"ellipse a/b={r:<4}    P={P:.6f}  I43={I43:.4f}  I32={I32:.4f}")
