"""
Étape 1H : Fit à 2 paramètres avec beta fixé.

Pour chaque p :
  1. beta_est = mediane(deficit) pour R >= 500
  2. resid(R) = deficit(R) - beta_est
  3. fit resid(R) = alpha * R^(-gamma) sur R in [100, 400]

Elimine la degenerescence du fit a 3 parametres.

Dépendances : numpy, matplotlib, scipy
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

try:
    # Import dynamique afin de permettre l'exécution avec les interpréteurs
    # dont l'analyseur statique ne connaît pas encore l'environnement Python.
    import importlib
    np = importlib.import_module("numpy")
except ImportError as exc:
    raise ImportError(
        "NumPy est requis. Installez-le avec : python -m pip install numpy"
    ) from exc
try:
    plt = importlib.import_module("matplotlib.pyplot")
except ImportError:
    class _MissingMatplotlib:
        def __getattr__(self, name):
            raise ImportError(
                "matplotlib est requis pour les graphiques. "
                "Installez-le avec : python -m pip install matplotlib"
            )

        def subplots(self, *args, **kwargs):
            raise ImportError(
                "matplotlib est requis pour les graphiques. "
                "Installez-le avec : python -m pip install matplotlib"
            )

        def savefig(self, *args, **kwargs):
            raise ImportError(
                "matplotlib est requis pour les graphiques. "
                "Installez-le avec : python -m pip install matplotlib"
            )

        def show(self, *args, **kwargs):
            return None

    plt = _MissingMatplotlib()
try:
    # Import dynamique pour éviter les faux positifs de l'analyseur statique
    # lorsque SciPy est installé dans un autre interpréteur Python.
    ConvexHull = importlib.import_module("scipy.spatial").ConvexHull
    binary_erosion = importlib.import_module("scipy.ndimage").binary_erosion
    linregress = importlib.import_module("scipy.stats").linregress
except ImportError as exc:
    raise ImportError(
        "SciPy est requis. Installez-le avec : python -m pip install scipy"
    ) from exc
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
# 2. STRATEGIE DE FIT EN 2 ETAPES
# ============================================================

def estimate_beta(R_arr, deficit, R_min_beta=500):
    """Mediane du deficit pour R >= R_min_beta."""
    R_arr = np.asarray(R_arr, dtype=float)
    deficit = np.asarray(deficit, dtype=float)
    mask = (R_arr >= R_min_beta) & np.isfinite(R_arr) & np.isfinite(deficit)
    if mask.sum() < 3:
        return np.nan
    return float(np.median(deficit[mask]))


def fit_alpha_gamma(R_arr, resid, R_min=100, R_max=400):
    """Fit resid = alpha * R^(-gamma) sur [R_min, R_max].

    Le residual brut deficit - beta peut devenir légèrement négatif au
    voisinage de la transition, alors que le modèle physique reste une loi
    de puissance sur la magnitude du residual. On essaie d'abord les points
    strictement positifs, puis la magnitude absolue pour éviter des NaN
    artificiels. 
    """
    R_arr = np.asarray(R_arr, dtype=float)
    resid = np.asarray(resid, dtype=float)

    mask = (R_arr >= R_min) & (R_arr <= R_max) & np.isfinite(R_arr) & np.isfinite(resid)
    if mask.sum() < 5:
        return np.nan, np.nan, np.nan

    R_fit = R_arr[mask]
    d_fit = resid[mask]

    for use_abs in (False, True):
        values = np.abs(d_fit) if use_abs else d_fit
        valid = values > 1e-15
        if valid.sum() < 5:
            continue

        logR = np.log(R_fit[valid])
        logD = np.log(values[valid])
        if not np.all(np.isfinite(logR)) or not np.all(np.isfinite(logD)):
            continue

        A = np.vstack([logR, np.ones_like(logR)]).T
        coef, _, _, _ = np.linalg.lstsq(A, logD, rcond=None)
        pred = A @ coef
        ss_res = np.sum((logD - pred) ** 2)
        ss_tot = np.sum((logD - logD.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return float(np.exp(coef[1])), float(-coef[0]), float(r2)

    return np.nan, np.nan, np.nan


# ============================================================
# 3. PROGRAMME PRINCIPAL
# ============================================================

def main():
    # Grille en p : fine autour de 2, grossiere ailleurs
    p_values = np.concatenate([
        np.arange(1.90, 2.00, 0.05),      # 1.90, 1.95
        np.arange(2.00, 2.50, 0.025),     # 2.000, 2.025, ..., 2.475
        np.arange(2.50, 3.00, 0.05),      # 2.50, 2.55, ..., 2.95
        [3.0, 4.0, 5.0, 7.0, 10.0],
    ])
    p_values = np.unique(np.round(p_values, 3))

    # Grille R logarithmique plus dense
    R_values = np.unique(
        np.round(np.geomspace(100, 800, 40)).astype(int)
    ).astype(float)

    print("=" * 70)
    print("FIT A 2 PARAMETRES AVEC BETA FIXE")
    print("=" * 70)
    print(f"p : {len(p_values)} valeurs de {p_values[0]:.2f} "
          f"a {p_values[-1]:.2f}")
    print(f"R : {len(R_values)} valeurs de {int(R_values[0])} "
          f"a {int(R_values[-1])}")
    print(f"Fenetre beta   : R >= 500")
    print(f"Fenetre fit    : R in [100, 400]")
    print()

    results = []
    curves = {}

    t0 = time.time()
    for i, p in enumerate(p_values):
        defs = np.empty(len(R_values))
        for j, R in enumerate(R_values):
            mask = superellipse_mask(R, p)
            Ph = hull_perimeter(mask)
            Pc = perimeter_continuous(R, p)
            defs[j] = abs(1.0 - Ph / Pc) if Pc > 0 else 0.0

        beta = estimate_beta(R_values, defs, R_min_beta=500)
        resid = defs - beta
        alpha, gamma, r2 = fit_alpha_gamma(R_values, resid,
                                            R_min=100, R_max=400)

        results.append((p, alpha, gamma, beta, r2))
        curves[p] = (R_values.copy(), defs.copy())

        print(f"  p = {p:5.3f}  |  "
              f"alpha = {alpha:>10.4e}  gamma = {gamma:>7.4f}  "
              f"beta = {beta:>10.4e}  R2 = {r2:.4f}  "
              f"({time.time() - t0:.1f}s)")

    # ============================================================
    # TABLEAU
    # ============================================================
    print("\n" + "=" * 70)
    print("RESULTATS")
    print("=" * 70)
    print(f"{'p':>7}  {'alpha':>12}  {'gamma':>9}  "
          f"{'beta':>12}  {'R2':>8}")
    print("-" * 60)
    for p, a, g, b, r2 in results:
        print(f"{p:>7.3f}  {a:>12.4e}  {g:>9.4f}  "
              f"{b:>12.4e}  {r2:>8.4f}")

    # ============================================================
    # TRACES
    # ============================================================
    p_arr = np.array([r[0] for r in results])
    alpha_arr = np.array([r[1] for r in results])
    gamma_arr = np.array([r[2] for r in results])
    beta_arr = np.array([r[3] for r in results])
    r2_arr = np.array([r[4] for r in results])

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # (1) beta(p)
    ax = axes[0, 0]
    ax.plot(p_arr, beta_arr, 'o-', color='b', markersize=5,
            linewidth=1.0)
    ax.axvline(2.0, color='r', linestyle='--', linewidth=1.2,
               label='p = 2')
    ax.set_xlabel('p')
    ax.set_ylabel('beta(p)')
    ax.set_title('Plateau asymptotique beta(p) [median R>=500]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (2) log-log beta vs (p-2)
    ax = axes[0, 1]
    mask = (p_arr > 2.0) & (beta_arr > 0)
    if mask.sum() > 3:
        dp = p_arr[mask] - 2.0
        logdp = np.log(dp)
        logbeta = np.log(beta_arr[mask])
        A = np.vstack([logdp, np.ones_like(logdp)]).T
        coef, _, _, _ = np.linalg.lstsq(A, logbeta, rcond=None)
        slope, intercept = coef
        ax.loglog(dp, beta_arr[mask], 'o', color='b', markersize=6,
                  label=f'beta (pente = {slope:.3f})')
        dp_fit = np.geomspace(dp.min(), dp.max(), 100)
        ax.loglog(dp_fit, np.exp(intercept) * dp_fit ** slope, '-',
                  color='b', linewidth=1.2,
                  label=f'fit : beta ~ (p-2)^{slope:.3f}')

        # Test d'une pente 1 (lineaire)
        ax.loglog(dp_fit, beta_arr[mask][0] * (dp_fit / dp[0]) ** 1.0,
                  '--', color='g', linewidth=1.0, label='pente 1')
        # Test d'une pente 2
        ax.loglog(dp_fit, beta_arr[mask][0] * (dp_fit / dp[0]) ** 2.0,
                  ':', color='r', linewidth=1.0, label='pente 2')
    ax.set_xlabel('p - 2')
    ax.set_ylabel('beta(p)')
    ax.set_title('Loi de puissance beta ~ (p-2)^s')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # (3) gamma(p)
    ax = axes[1, 0]
    ax.plot(p_arr, gamma_arr, 's-', color='g', markersize=5,
            linewidth=1.0)
    ax.axvline(2.0, color='r', linestyle='--', linewidth=1.2)
    ax.axhline(4/3, color='k', linestyle=':', linewidth=1.0,
               label='gamma = 4/3')
    ax.set_xlabel('p')
    ax.set_ylabel('gamma(p)')
    ax.set_title('Exposant gamma(p) — fit a 2 parametres')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (4) R^2 du fit
    ax = axes[1, 1]
    ax.plot(p_arr, r2_arr, 'd-', color='m', markersize=5,
            linewidth=1.0)
    ax.axvline(2.0, color='r', linestyle='--', linewidth=1.2)
    ax.axhline(0.99, color='k', linestyle=':', linewidth=1.0,
               label='R2 = 0.99')
    ax.set_xlabel('p')
    ax.set_ylabel('R^2 du fit log-log')
    ax.set_title('Qualite du fit (log-log)')
    ax.set_ylim(0.9, 1.01)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('transition_p2_v2.png', dpi=150)
    try:
        plt.show()
    except Exception:
        pass

    # ============================================================
    # ANALYSE
    # ============================================================
    print("\n" + "=" * 70)
    print("ANALYSE DE LA TRANSITION")
    print("=" * 70)

    # Localisation du saut
    beta_valid = np.isfinite(beta_arr)
    if beta_valid.sum() > 1:
        jumps = np.abs(np.diff(beta_arr[beta_valid]))
        if jumps.size > 0:
            i_jump = np.argmax(jumps)
            i_jump_global = np.where(beta_valid)[0][i_jump]
            print(f"Plus grand saut de beta : entre p = {p_arr[i_jump_global]:.3f} "
                  f"et p = {p_arr[i_jump_global + 1]:.3f}")
            print(f"  beta = {beta_arr[i_jump_global]:.4e} -> "
                  f"{beta_arr[i_jump_global + 1]:.4e}")

    # Fit beta ~ (p-2)^s
    mask = (p_arr > 2.02) & np.isfinite(beta_arr) & (beta_arr > 0)
    if mask.sum() > 3:
        dp = p_arr[mask] - 2.0
        beta_valid_data = beta_arr[mask]
        s = np.polyfit(np.log(dp), np.log(beta_valid_data), 1)
        print(f"\nFit log-log beta ~ (p-2)^s")
        print(f"  s = {s[0]:.4f}")
        for s_test in [1, 2, 3, 4]:
            log_dp = np.log(dp)
            log_beta = np.log(beta_valid_data)
            pred = s_test * log_dp + np.mean(log_beta - s_test * log_dp)
            ss_res = np.sum((log_beta - pred) ** 2)
            ss_tot = np.sum((log_beta - log_beta.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
            print(f"  test s = {s_test} : R2 = {r2:.5f}")

    # gamma moyen avant / apres
    mask_le2 = (p_arr <= 2.0) & np.isfinite(gamma_arr)
    mask_gt22 = (p_arr > 2.2) & np.isfinite(gamma_arr)
    if mask_le2.sum() > 0:
        print(f"\ngamma moyen p <= 2.00 : "
              f"{np.nanmean(gamma_arr[mask_le2]):.4f}")
    if mask_gt22.sum() > 0:
        print(f"gamma moyen p >= 2.20 : "
              f"{np.nanmean(gamma_arr[mask_gt22]):.4f}")

    # Verifier que beta ~ 0 pour p <= 2
    beta_petit = beta_arr[(p_arr <= 2.0) & np.isfinite(beta_arr)]
    if beta_petit.size > 0:
        print(f"\nbeta pour p <= 2 : max = {beta_petit.max():.3e}")
    else:
        print("\nbeta pour p <= 2 : aucune valeur finie")

    print("""
INTERPRETATION
--------------
1. Si gamma(p) reste stable autour de 1.27 pour tout p, la
   structure du deficit est universelle : un plateau beta(p)
   plus une loi de puissance de meme exposant.

2. Si gamma(p) augmente avec p, le plateau est atteint plus
   rapidement et la partie transitoire est plus raide.

3. L'exposant critique s de beta ~ (p-2)^s est determine par
   comparaison des R^2 pour s = 1, 2, 3, 4. Si le meilleur
   fit est clairement s = 2, c'est une transition du second
   ordre au sens de Landau.

4. La valeur de R^2 indique si le modele a 2 parametres est
   adequat. Si R^2 > 0.99 partout, le fit est fiable.
""")

    print(f"\nTemps total : {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()