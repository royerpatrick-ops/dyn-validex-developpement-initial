"""
Modèle local : enveloppe convexe d'un réseau aléatoire (mesure de Haar) sous la parabole y = -x²/2.

Après l'application affine qui conserve les aires (x' = κ^{1/3} x_t, y' = κ^{-1/3} y_n),
le bord de RK devient localement la parabole y = -x²/2 et Z² un réseau unimodulaire.
Dans ce repère :
  - densité de sommets par unité de longueur      -> doit valoir c0/2π ≈ 0,5496
  - densité de déficit D = < g - g'^2/2 >          -> doit valoir C∞
où g(x) = -x²/2 - h(x) est l'écart entre la parabole et l'enveloppe supérieure h.
(Déficit d'un contour intérieur à distance δ(s) : ∫(κδ - δ'^2/2) ds au second ordre.)

Échantillonnage de Haar sur SL2(R)/SL2(Z) : point z = x + i y du domaine fondamental
(densité dx dy / y², masse totale π/3), rotation θ uniforme, translation uniforme.
Réseau : v1 = (1/√y)(1, 0), v2 = (1/√y)(x, y), tournés de θ ; |v1| = 1/√y est le plus court.
Stratification en y : tranche basse [√3/2, 1], puis tranches dyadiques [2^k, 2^{k+1}].
"""
import math, sys, time
import numpy as np  # pyright: ignore[reportMissingImports]
from scipy.spatial import ConvexHull  # pyright: ignore[reportMissingImports]

GL_X = np.array([-math.sqrt(3 / 5), 0.0, math.sqrt(3 / 5)])
GL_W = np.array([5 / 9, 8 / 9, 5 / 9])


def poids_tranche(a, b):
    """Mesure de Haar (dx dy / y²) de F ∩ {a ≤ y < b}."""
    xs = np.linspace(-0.5, 0.5, 20001)
    y0 = np.sqrt(1 - xs ** 2)
    lo = np.maximum(a, y0)
    f = np.where(lo < b, 1 / lo - 1 / b, 0.0)
    return float(np.trapezoid(f, xs)) if hasattr(np, "trapezoid") else float(np.trapz(f, xs))


def tirer_z(a, b, rng):
    """Tire (x, y) dans F ∩ {a ≤ y < b} avec densité ∝ 1/y²."""
    while True:
        x = rng.uniform(-0.5, 0.5)
        y0 = max(a, math.sqrt(1 - x * x))
        if y0 >= b:
            continue
        # y de densité ∝ 1/y² sur [y0, b) ; on accepte selon la masse de la colonne
        U = rng.uniform()
        y = 1 / (1 / y0 - U * (1 / y0 - 1 / b))
        masse = 1 / y0 - 1 / b
        masse_max = 1 / max(a, math.sqrt(3) / 2) - 1 / b
        if rng.uniform() * masse_max <= masse:
            return x, y


def enveloppe(x, y, th, u, v, W):
    """Points extrêmes de chaque droite du réseau (parallèle à v1) sous la parabole,
    dans la bande |X| ≤ W ; renvoie la chaîne supérieure de l'enveloppe (X croissants)."""
    s = 1 / math.sqrt(y)
    c, si = math.cos(th), math.sin(th)
    rot = lambda p: np.array([c * p[0] - si * p[1], si * p[0] + c * p[1]])
    v1 = rot(np.array([s, 0.0])); v2 = rot(np.array([x * s, y * s]))
    o = u * v1 + v * v2
    n = np.array([-v1[1], v1[0]]) / s                       # normale unitaire aux droites
    pas = abs(n @ v2)                                       # = √y, écart entre droites
    # boîte qui contient la région utile
    coins = np.array([[-W, 0.5], [W, 0.5], [-W, -W * W / 2 - 1], [W, -W * W / 2 - 1]])
    proj = coins @ n - n @ o
    sg = np.sign(n @ v2)
    bmin = math.floor(proj.min() / pas) - 1; bmax = math.ceil(proj.max() / pas) + 1
    b = np.arange(bmin, bmax + 1) * sg
    cx = o[0] + b * v2[0]; cy = o[1] + b * v2[1]
    # p(a) = c + a v1 ; contrainte cy + a v1y + (cx + a v1x)²/2 ≤ 0
    A = v1[0] ** 2 / 2; B = v1[1] + cx * v1[0]; C = cy + cx ** 2 / 2
    # racines de A a² + B a + C = 0 par la formule stable (sans seuil sur A)
    disc = B * B - 4 * A * C
    ok = disc >= 0
    r = np.sqrt(np.where(ok, disc, 0.0))
    q = -0.5 * (B + np.where(B >= 0, r, -r))
    with np.errstate(divide="ignore", invalid="ignore"):
        r1 = np.where(A > 0, q / A, np.where(B > 0, -np.inf, np.inf))
        r2 = np.where(q != 0, C / q, np.where(B > 0, -np.inf, np.inf))
    a1 = np.minimum(r1, r2); a2 = np.maximum(r1, r2)
    if A == 0:  # v1 exactement vertical : contrainte linéaire B a + C <= 0
        a1 = np.where(B > 0, -np.inf, -C / B); a2 = np.where(B > 0, -C / B, np.inf)
    # contrainte |X| ≤ W
    if abs(v1[0]) > 1e-15:
        t1 = (-W - cx) / v1[0]; t2 = (W - cx) / v1[0]
        a1 = np.maximum(a1, np.minimum(t1, t2)); a2 = np.minimum(a2, np.maximum(t1, t2))
    # contrainte de profondeur Y ≥ -W²/2 - 1 (borne aussi les droites quasi verticales)
    if abs(v1[1]) > 0:
        t1 = (-W * W / 2 - 1 - cy) / v1[1]; t2 = (0.5 - cy) / v1[1]
        a1 = np.maximum(a1, np.minimum(t1, t2)); a2 = np.minimum(a2, np.maximum(t1, t2))
    lo = np.ceil(a1); hi = np.floor(a2)
    dehors = lambda aa: cy + aa * v1[1] + (cx + aa * v1[0]) ** 2 / 2 > 0
    for _ in range(3):
        lo = np.where(ok & (lo <= hi) & dehors(lo), lo + 1, lo)
        hi = np.where(ok & (lo <= hi) & dehors(hi), hi - 1, hi)
    ok &= lo <= hi
    pts = []
    for aa in (lo[ok], hi[ok]):
        pts.append(np.column_stack([cx[ok] + aa * v1[0], cy[ok] + aa * v1[1]]))
    P = np.vstack(pts)
    # sécurité : retirer d'éventuels points au-dessus de la parabole (arrondis)
    P = P[P[:, 1] <= -P[:, 0] ** 2 / 2 + 1e-12]
    if len(P) < 3:
        return None
    h = ConvexHull(P)
    V = P[h.vertices]                                       # sens trigonométrique
    i0 = np.argmax(V[:, 0]); V = np.roll(V, -i0, axis=0)    # part du point le plus à droite
    i1 = np.argmin(V[:, 0])
    haut = V[: i1 + 1][::-1]                                # chaîne supérieure, X croissants
    return haut


def mesures(haut, Wc):
    """Sur |X| ≤ Wc : ∫(g - g'²/2) dX, nombre de sommets, et arêtes entièrement incluses."""
    X, Y = haut[:, 0], haut[:, 1]
    D = 0.0
    for i in range(len(X) - 1):
        x1, x2 = max(X[i], -Wc), min(X[i + 1], Wc)
        if x2 <= x1:
            continue
        pente = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])
        xm, dx = (x1 + x2) / 2, (x2 - x1) / 2
        t = xm + dx * GL_X
        g = -t ** 2 / 2 - (Y[i] + pente * (t - X[i]))
        gp = -t - pente
        D += dx * np.sum(GL_W * (g - gp ** 2 / 2))
    inside = (X >= -Wc) & (X <= Wc)
    nsom = int(inside.sum())
    j = np.where(inside[:-1] & inside[1:])[0]
    w = X[j + 1] - X[j]
    contrib = []
    for i in j:
        pente = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])
        xm, dx = (X[i] + X[i + 1]) / 2, (X[i + 1] - X[i]) / 2
        t = xm + dx * GL_X
        g = -t ** 2 / 2 - (Y[i] + pente * (t - X[i])); gp = -t - pente
        contrib.append(dx * np.sum(GL_W * (g - gp ** 2 / 2)))
    return D, nsom, w, np.array(contrib)


def main(n_par_tranche=20000, kmax=24, graine=1):
    """Mesure ponctuelle : fenêtre centrale |X| ≤ 1/2 (longueur 1), pour que la tranche
    en y décrive bien le réseau local en X = 0. Fenêtre de calcul W adaptée à la tranche."""
    rng = np.random.default_rng(graine)
    tranches = [(math.sqrt(3) / 2, 1.0)] + [(2.0 ** k, 2.0 ** (k + 1)) for k in range(kmax)]
    total = math.pi / 3
    lignes = []; AW, AC, AP = [], [], []
    t0 = time.time()
    for (a, b) in tranches:
        p = poids_tranche(a, b) / total
        W = 8 + 2 * max(6 * b ** 0.25, math.sqrt(b))
        Ds, Ns = np.empty(n_par_tranche), np.empty(n_par_tranche)
        for k in range(n_par_tranche):
            x, y = tirer_z(a, b, rng)
            th = rng.uniform(0, 2 * math.pi); u, v = rng.uniform(0, 1, 2)
            haut = enveloppe(x, y, th, u, v, W)
            _, ns, _, _ = mesures(haut, 0.5)
            Ns[k] = ns; Dk = 0.0
            # arêtes dont le milieu est dans [-1/2, 1/2] : longueur et contribution totale
            X, Y = haut[:, 0], haut[:, 1]; mid = (X[:-1] + X[1:]) / 2
            for i in np.where(np.abs(mid) <= 0.5)[0]:
                pente = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])
                xm, dx = mid[i], (X[i + 1] - X[i]) / 2
                t = xm + dx * GL_X
                g = -t ** 2 / 2 - (Y[i] + pente * (t - X[i])); gp = -t - pente
                ce = dx * np.sum(GL_W * (g - gp ** 2 / 2)); Dk += ce
                AW.append(2 * dx); AC.append(ce); AP.append(p / n_par_tranche)
            Ds[k] = Dk
        se = Ds.std() / math.sqrt(n_par_tranche)
        lignes.append((a, b, p, Ds.mean(), se, Ns.mean(), Ns.std() / math.sqrt(n_par_tranche)))
        print(f"  y ∈ [{a:9.3f}, {b:9.0f})  poids {p:.3e}  D = {Ds.mean():9.4f} ± {se:.4f}"
              f"   p·D = {p*Ds.mean():.5f} ± {p*se:.5f}   sommets = {Ns.mean():.3f}   ({time.time()-t0:.0f} s)", flush=True)
    L = np.array(lignes)
    np.save("modele_local.npy", {"tranches": L, "w": np.array(AW), "c": np.array(AC), "p": np.array(AP)},
            allow_pickle=True)
    pD = L[:, 2] * L[:, 3]; pN = L[:, 2] * L[:, 5]
    print(f"\nSomme sur les tranches calculées : D = {pD.sum():.4f} ± {np.sqrt(np.sum((L[:,2]*L[:,4])**2)):.4f}"
          f" ; sommets par unité de longueur = {pN.sum():.4f} ± {np.sqrt(np.sum((L[:,2]*L[:,6])**2)):.4f}")
    return L


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    main(n)
