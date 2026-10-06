"""
Étape 1I : Exposant local gamma_eff(R) par regression glissante.

Pour chaque p :
  - Calcul du deficit(R) sur une grille fine en log R
  - Estimation locale de gamma_eff(R) par regression lineaire
    sur une fenetre glissante en log-log
  - Extraction de gamma_0 (valeur a petit R) et de R_c (croisement)

Dépendances : numpy, matplotlib, scipy
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import numpy as np  # type: ignore[import-not-found]
try:
    import matplotlib.pyplot as plt  # type: ignore[import-not-found]
except Exception:  # pragma: no cover - matplotlib is optional for CI/static checks
    plt = None
try:
    from scipy.spatial import ConvexHull  # type: ignore[import-not-found, reportMissingImports]
except Exception:  # pragma: no cover - scipy is optional in some environments
    ConvexHull = None
try:
    from scipy.ndimage import binary_erosion  # type: ignore[import-not-found, reportMissingImports]
except Exception:  # pragma: no cover - scipy is optional in some environments
    binary_erosion = None
import time

_integrate = getattr(np, "trapezoid", None) or np.trapz

THETA0 = 0.3


# ============================================================
# 1. GEOMETRIE
# ============================================================

def superellipse_mask(R, p, theta0=THETA0):
    L = int(np.ceil(R)) + 1
    x = np.arange(-L, L + 1)
    X, Y = np.meshgrid(x, x)
    c, s = np.cos(theta0), np.sin(theta0)
    Xr = c * X + s * Y
    Yr = -s * X + c * Y
    return (np.abs(Xr / R) ** p + np.abs(Yr / R) ** p) <= 1.0


def hull_perimeter(mask):
    eroded = binary_erosion(mask, iterations=1)
    border = mask & ~eroded
    ys, xs = np.where(border)
    if len(xs) < 3:
        return 0.0
    pts = np.column_stack([xs, ys]).astype(float)
    try:
        hull = ConvexHull(pts)
    except Exception:
        return 0.0
    v = pts[hull.vertices]
    d = np.diff(np.vstack([v, v[:1]]), axis=0)
    return float(np.sqrt((d ** 2).sum(axis=1)).sum())


def perimeter_continuous(R, p, N_theta=20000):
    if abs(p - 1.0) < 1e-9:
        return 4 * np.sqrt(2) * R
    if p > 1e8:
        return 8 * R
    q = p / (p - 1)
    theta = np.linspace(0, np.pi, N_theta)
    integrand = (np.abs(np.cos(theta)) ** q
                 + np.abs(np.sin(theta)) ** q) ** (1.0 / q)
    return 2 * R * _integrate(integrand, theta)


# ============================================================
# 2. EXPOSANT LOCAL PAR REGRESSION GLISSANTE
# ============================================================

def gamma_local(R, deficit, window):
    """
    gamma_eff(R_i) = pente locale de log(deficit) vs log(R),
    estimee par regression lineaire sur une fenetre [i-w, i+w].

    Retourne un tableau de meme taille que R, avec NaN aux bords.
    """
    logR = np.log(R)
    logD = np.log(deficit)
    n = len(R)
    g = np.full(n, np.nan)
    for i in range(window, n - window):
        lo = i - window
        hi = i + window + 1
        x = logR[lo:hi]
        y = logD[lo:hi]
        # Regression lineaire
        A = np.vstack([x, np.ones_like(x)]).T
        coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        g[i] = -coef[0]
    return g


def extract_Rc(R, gamma_eff, gamma_0, threshold=0.5):
    """
    R_c : premier R ou gamma_eff(R) < threshold * gamma_0.
    """
    target = threshold * gamma_0
    for i in range(len(R)):
        if np.isfinite(gamma_eff[i]) and gamma_eff[i] < target:
            return R[i]
    return np.nan


# ============================================================
# 3. PROGRAMME PRINCIPAL
# ============================================================

def main():
    # p selectionnes pour couvrir les regimes
    p_values = [1.5, 1.9, 2.0, 2.05, 2.1, 2.15,
                2.2, 2.3, 2.5, 3.0, 5.0, 10.0]

    # Grille R dense en log
    R_values = np.unique(
        np.round(np.geomspace(80, 1000, 60)).astype(int)
    ).astype(float)

    # Taille de la fenetre glissante (en nombre de points)
    window = 3

    print("=" * 70)
    print("EXPOSANT LOCAL gamma_eff(R)")
    print("=" * 70)
    print(f"p : {len(p_values)} valeurs : {p_values}")
    print(f"R : {len(R_values)} valeurs de {int(R_values[0])} "
          f"a {int(R_values[-1])}")
    print(f"Fenetre glissante : +/- {window} points")
    print()

    results = {}
    gamma_curves = {}

    t0 = time.time()
    for p in p_values:
        # Calcul du deficit
        defs = np.empty(len(R_values))
        for j, R in enumerate(R_values):
            mask = superellipse_mask(R, p)
            Ph = hull_perimeter(mask)
            Pc = perimeter_continuous(R, p)
            defs[j] = abs(1.0 - Ph / Pc) if Pc > 0 else 0.0

        # Exposant local
        g_eff = gamma_local(R_values, defs, window)

        # gamma_0 : mediane de gamma_eff sur les premiers points valides
        valid = np.isfinite(g_eff)
        if valid.sum() > 5:
            # Prendre les 20% premiers points valides
            idx_valid = np.where(valid)[0]
            n_first = max(3, len(idx_valid) // 5)
            gamma_0 = float(np.median(g_eff[idx_valid[:n_first]]))
        else:
            gamma_0 = np.nan

        # R_c : croisement a gamma_0 / 2
        R_c = extract_Rc(R_values, g_eff, gamma_0, threshold=0.5)

        # gamma_final : mediane sur les derniers points valides
        if valid.sum() > 5:
            n_last = max(3, len(idx_valid) // 5)
            gamma_final = float(np.median(g_eff[idx_valid[-n_last:]]))
        else:
            gamma_final = np.nan

        results[p] = {
            'gamma_0': gamma_0,
            'gamma_final': gamma_final,
            'R_c': R_c,
        }
        gamma_curves[p] = (R_values.copy(), defs.copy(), g_eff.copy())

        print(f"  p = {p:5.2f}  |  "
              f"gamma_0 = {gamma_0:>6.3f}  "
              f"gamma_final = {gamma_final:>6.3f}  "
              f"R_c = {R_c if np.isfinite(R_c) else '---':>6}  "
              f"({time.time() - t0:.1f}s)")

    # ============================================================
    # TABLEAU RECAPITULATIF
    # ============================================================
    print("\n" + "=" * 70)
    print("RESULTATS")
    print("=" * 70)
    print(f"{'p':>6}  {'gamma_0':>9}  {'gamma_final':>12}  {'R_c':>8}")
    print("-" * 42)
    for p in p_values:
        r = results[p]
        R_c_str = f"{r['R_c']:.0f}" if np.isfinite(r['R_c']) else "---"
        print(f"{p:>6.2f}  {r['gamma_0']:>9.4f}  "
              f"{r['gamma_final']:>12.4f}  {R_c_str:>8}")

    # ============================================================
    # TRACES
    # ============================================================
    if plt is not None:
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))

        # (1) gamma_eff(R) pour tous les p
        ax = axes[0, 0]
        cmap = plt.cm.coolwarm
        for p in p_values:
            R_arr, _, g_eff = gamma_curves[p]
            color = cmap((p - p_values[0]) / (p_values[-1] - p_values[0]))
            ax.semilogx(R_arr, g_eff, '-', color=color,
                        linewidth=1.4, label=f'p = {p}')
        ax.axhline(4/3, color='k', linestyle=':', linewidth=1.0,
                   label='4/3')
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.5)
        ax.set_xlabel('R')
        ax.set_ylabel('gamma_eff(R)')
        ax.set_title('Exposant local — tous les p')
        ax.legend(fontsize=7, ncol=2, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.5, 3.5)

        # (2) Zoom sur les p de transition
        ax = axes[0, 1]
        p_zoom = [2.0, 2.05, 2.1, 2.15, 2.2, 2.3]
        colors_zoom = plt.cm.plasma(np.linspace(0, 0.9, len(p_zoom)))
        for k, p in enumerate(p_zoom):
            R_arr, _, g_eff = gamma_curves[p]
            ax.semilogx(R_arr, g_eff, '-', color=colors_zoom[k],
                        linewidth=1.8, label=f'p = {p}')
        ax.axhline(4/3, color='k', linestyle=':', linewidth=1.0)
        ax.axhline(0, color='gray', linewidth=0.5)
        ax.set_xlabel('R')
        ax.set_ylabel('gamma_eff(R)')
        ax.set_title('Zoom : zone de transition p ~ 2')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.5, 3.5)

        # (3) Deficit(R) pour tous les p
        ax = axes[1, 0]
        for p in p_values:
            R_arr, defs, _ = gamma_curves[p]
            color = cmap((p - p_values[0]) / (p_values[-1] - p_values[0]))
            ax.loglog(R_arr, defs, '-', color=color, linewidth=1.2)
        sm = plt.cm.ScalarMappable(cmap=cmap,
            norm=plt.Normalize(vmin=p_values[0], vmax=p_values[-1]))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('p')
        ax.set_xlabel('R')
        ax.set_ylabel('Deficit relatif')
        ax.set_title('Deficit(R) — vue brute')
        ax.grid(True, alpha=0.3, which='both')

        # (4) gamma_0 et R_c vs p
        ax = axes[1, 1]
        p_arr = np.array(p_values)
        g0_arr = np.array([results[p]['gamma_0'] for p in p_values])
        Rc_arr = np.array([results[p]['R_c'] for p in p_values])
        ax.plot(p_arr, g0_arr, 'o-', color='b', markersize=6,
                linewidth=1.2, label='gamma_0(p)')
        ax.axhline(4/3, color='k', linestyle=':', linewidth=1.0,
                   label='4/3')
        ax.set_xlabel('p')
        ax.set_ylabel('gamma_0', color='b')
        ax.tick_params(axis='y', labelcolor='b')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=9)
        ax2 = ax.twinx()
        valid_Rc = np.isfinite(Rc_arr)
        if valid_Rc.any():
            ax2.semilogy(p_arr[valid_Rc], Rc_arr[valid_Rc],
                         's--', color='r', markersize=6,
                         linewidth=1.2, label='R_c(p)')
            ax2.set_ylabel('R_c', color='r')
            ax2.tick_params(axis='y', labelcolor='r')
            ax2.legend(loc='upper right', fontsize=9)

        plt.tight_layout()
        plt.savefig('gamma_local.png', dpi=150)
        plt.show()
    else:
        print("Matplotlib non disponible : aucun graphe produit.")

    # ============================================================
    # ANALYSE
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSE")
    print("=" * 70)

    # gamma_0 pour p <= 2 vs p > 2
    mask_le2 = p_arr <= 2.0
    mask_gt2 = p_arr > 2.0
    if mask_le2.any():
        print(f"gamma_0 moyen pour p <= 2.00 : "
              f"{np.nanmean(g0_arr[mask_le2]):.4f}")
    if mask_gt2.any():
        print(f"gamma_0 moyen pour p >  2.00 : "
              f"{np.nanmean(g0_arr[mask_gt2]):.4f}")

    # Verifier que R_c est bien defini
    if valid_Rc.any():
        Rc_valid = Rc_arr[valid_Rc]
        print(f"\nR_c : min = {Rc_valid.min():.0f}  "
              f"max = {Rc_valid.max():.0f}  "
              f"median = {np.median(Rc_valid):.0f}")

    # Detection de la decroissance de gamma_eff pour p > 2
    print("\nDecroissance de gamma_eff (ratio gamma_final/gamma_0) :")
    for p in p_values:
        r = results[p]
        if np.isfinite(r['gamma_0']) and r['gamma_0'] > 0:
            ratio = r['gamma_final'] / r['gamma_0']
            label = "plateau" if ratio < 0.5 else "puissance"
            print(f"  p = {p:>5.2f} : ratio = {ratio:>6.3f}  -> {label}")

    print("""
INTERPRETATION
--------------
1. Pour p <= 2 (courbure stricte), gamma_eff(R) doit rester
   approximativement constant autour de 1.27-1.33 : le deficit
   suit une loi de puissance pure.

2. Pour p > 2 (points plats), gamma_eff(R) doit DECROITRE de
   gamma_0 vers 0 : c'est la signature d'un plateau asymptotique.
   Le croisement a gamma_0/2 definit R_c(p).

3. Si R_c(p) diminue quand p augmente, la transition vers le
   plateau est plus precoce. Si R_c(p) augmente, c'est l'inverse.

4. Si gamma_0 varie peu avec p (par exemple 1.27-1.40), alors la
   partie loi-de-puissance du deficit est UNIVERSLLE, et seul
   beta(p) distingue les regimes. C'est la conclusion la plus
   propre pour la conjecture.
""")

    print(f"\nTemps total : {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()