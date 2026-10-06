"""
Étape 1B : Test de fenêtre glissante de l'exposant γ.

Pour chaque R0, on ajuste |π_eff(R) - π| ~ α R^{-γ} sur une fenêtre
[R0, R0+L] (linéaire) ou [R0, c·R0] (logarithmique), et on trace
γ_eff(R0). Révèle si γ converge lentement vers la valeur asymptotique
(3/2 pour l'aire, 4/3 pour le périmètre) ou s'il est multifractal.

Dépendances : numpy, matplotlib
"""

import importlib

try:
    np = importlib.import_module("numpy")
except ImportError as exc:
    raise ImportError(
        "NumPy est requis. Installez-le avec : python -m pip install numpy"
    ) from exc
try:
    plt = importlib.import_module("matplotlib.pyplot")
except ImportError as exc:
    raise ImportError(
        "Matplotlib est requis. Installez-le avec : python -m pip install matplotlib"
    ) from exc
import math
import time

# ============================================================
# 1. CALCULS GÉOMÉTRIQUES
# ============================================================

def lattice_count(R):
    """Nombre de points entiers (x,y) avec x²+y² ≤ R²."""
    c = 0
    for x in range(-R, R + 1):
        c += 2 * math.isqrt(R * R - x * x) + 1
    return c

def hull_perimeter(R):
    """
    Périmètre de l'enveloppe convexe des points entiers du disque
    de rayon R, par monotone chain sur y_top(x) = ⌊√(R²−x²)⌋.
    Complexité O(R).
    """
    xs = np.arange(-R, R + 1)
    y_top = np.array([math.isqrt(R * R - int(x) * int(x)) for x in xs],
                     dtype=float)
    upper = []
    for x, y in zip(xs, y_top):
        while len(upper) >= 2:
            x1, y1 = upper[-2]
            x2, y2 = upper[-1]
            cross = (x2 - x1) * (y - y2) - (y2 - y1) * (x - x2)
            if cross >= 0:
                upper.pop()
            else:
                break
        upper.append((int(x), int(y)))
    lower_rev = [(x, -y) for (x, y) in reversed(upper)]
    hull = np.array(upper + lower_rev[1:-1], dtype=float)
    diffs = np.diff(np.vstack([hull, hull[:1]]), axis=0)
    return float(np.sum(np.sqrt((diffs ** 2).sum(axis=1))))

# ============================================================
# 2. AJUSTEMENT LOCAL DE γ
# ============================================================

def fit_gamma_local(R, dev, i0, i1):
    """Ajuste γ par régression log-log linéaire sur l'intervalle [i0, i1]."""
    mask = slice(i0, i1)
    x = R[mask]
    y = dev[mask]
    valid = y > 0
    if np.sum(valid) < 15:
        return np.nan, np.nan
    x = x[valid]
    y = y[valid]
    logx = np.log(x)
    logy = np.log(y)
    A = np.vstack([logx, np.ones_like(logx)]).T
    coef, _, _, _ = np.linalg.lstsq(A, logy, rcond=None)
    slope, intercept = coef
    return np.exp(intercept), -slope

def sliding_gamma_linear(R, dev, L, min_pts=15):
    """γ_eff(R0) sur fenêtre [R0, R0+L]."""
    n = len(R)
    gamma = np.full(n, np.nan)
    for i in range(n):
        R1 = R[i] + L
        j = np.searchsorted(R, R1)
        if j - i < min_pts:
            continue
        _, g = fit_gamma_local(R, dev, i, j)
        gamma[i] = g
    return gamma

def sliding_gamma_log(R, dev, c, min_pts=15):
    """γ_eff(R0) sur fenêtre [R0, c·R0]."""
    n = len(R)
    gamma = np.full(n, np.nan)
    for i in range(n):
        R1 = c * R[i]
        j = np.searchsorted(R, R1)
        if j - i < min_pts:
            continue
        _, g = fit_gamma_local(R, dev, i, j)
        gamma[i] = g
    return gamma

# ============================================================
# 3. PROGRAMME PRINCIPAL
# ============================================================

def main():
    R_max = 1200
    R_values = np.arange(5, R_max + 1)
    n = len(R_values)
    R = R_values.astype(float)

    print("=" * 62)
    print("TEST DE FENÊTRE GLISSANTE DE L'EXPOSANT γ")
    print("=" * 62)

    # --- Calcul des observables ---
    areas = np.empty(n)
    perimeters = np.empty(n)
    t0 = time.time()
    for i, r in enumerate(R_values):
        areas[i] = lattice_count(int(r))
        perimeters[i] = hull_perimeter(int(r))
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{n}   {time.time()-t0:.1f}s")

    pi_area = areas / R ** 2
    pi_perim = perimeters / (2 * R)

    dev_area = np.abs(pi_area - np.pi)
    dev_perim = np.abs(pi_perim - np.pi)

    # ============================================================
    # FENÊTRES GLISSANTES
    # ============================================================
    print("\nCalcul des fenêtres glissantes...")

    Ls = [50, 100, 200, 400]
    cs = [1.2, 1.5, 2.0]

    gamma_area_lin = {L: sliding_gamma_linear(R, dev_area, L) for L in Ls}
    gamma_perim_lin = {L: sliding_gamma_linear(R, dev_perim, L) for L in Ls}
    gamma_area_log = {c: sliding_gamma_log(R, dev_area, c) for c in cs}
    gamma_perim_log = {c: sliding_gamma_log(R, dev_perim, c) for c in cs}

    # --- Valeurs à R = 400 et R = 800 ---
    print("\n" + "=" * 62)
    print("VALEURS DE γ_eff EN FONCTION DE L'ÉCHELLE")
    print("=" * 62)

    for R_probe in [200, 400, 800]:
        i = np.argmin(np.abs(R - R_probe))
        print(f"\nÀ R = {int(R[i])} :")
        print(f"  Aire      : γ_log(1.5) = {gamma_area_log[1.5][i]:.4f}, "
              f"γ_lin(200) = {gamma_area_lin[200][i]:.4f}")
        print(f"  Périmètre : γ_log(1.5) = {gamma_perim_log[1.5][i]:.4f}, "
              f"γ_lin(200) = {gamma_perim_lin[200][i]:.4f}")

    # ============================================================
    # TRACÉS
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (1) Aire, fenêtres linéaires
    ax = axes[0, 0]
    for L in Ls:
        ax.plot(R, gamma_area_lin[L], '-', linewidth=1.2,
                label=f'L = {L}')
    ax.axhline(1.5, color='k', linestyle='--', linewidth=1.5,
               label='γ = 3/2 (asymptotique)')
    ax.set_xlabel('R₀')
    ax.set_ylabel('γ_eff(R₀)')
    ax.set_title("Aire — fenêtres linéaires [R₀, R₀+L]")
    ax.set_ylim(0.8, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (2) Aire, fenêtres logarithmiques
    ax = axes[0, 1]
    for c in cs:
        ax.plot(R, gamma_area_log[c], '-', linewidth=1.2,
                label=f'c = {c}')
    ax.axhline(1.5, color='k', linestyle='--', linewidth=1.5,
               label='γ = 3/2 (asymptotique)')
    ax.set_xlabel('R₀')
    ax.set_ylabel('γ_eff(R₀)')
    ax.set_title("Aire — fenêtres log [R₀, c·R₀]")
    ax.set_ylim(0.8, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (3) Périmètre, fenêtres linéaires
    ax = axes[1, 0]
    for L in Ls:
        ax.plot(R, gamma_perim_lin[L], '-', linewidth=1.2,
                label=f'L = {L}')
    ax.axhline(4/3, color='k', linestyle='--', linewidth=1.5,
               label='γ = 4/3 (asymptotique)')
    ax.set_xlabel('R₀')
    ax.set_ylabel('γ_eff(R₀)')
    ax.set_title("Périmètre — fenêtres linéaires [R₀, R₀+L]")
    ax.set_ylim(0.8, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (4) Périmètre, fenêtres logarithmiques
    ax = axes[1, 1]
    for c in cs:
        ax.plot(R, gamma_perim_log[c], '-', linewidth=1.2,
                label=f'c = {c}')
    ax.axhline(4/3, color='k', linestyle='--', linewidth=1.5,
               label='γ = 4/3 (asymptotique)')
    ax.set_xlabel('R₀')
    ax.set_ylabel('γ_eff(R₀)')
    ax.set_title("Périmètre — fenêtres log [R₀, c·R₀]")
    ax.set_ylim(0.8, 2.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('sliding_window_gamma.png', dpi=150)
    plt.show()

    # ============================================================
    # ANALYSE DE CONVERGENCE
    # ============================================================
    print("\n" + "=" * 62)
    print("ANALYSE DE CONVERGENCE DE γ_eff")
    print("=" * 62)

    for label, gammas, target in [
        ("Aire (log, c=1.5)", gamma_area_log[1.5], 1.5),
        ("Périmètre (log, c=1.5)", gamma_perim_log[1.5], 4/3),
        ("Aire (lin, L=200)", gamma_area_lin[200], 1.5),
        ("Périmètre (lin, L=200)", gamma_perim_lin[200], 4/3),
    ]:
        g = gammas[~np.isnan(gammas)]
        if len(g) < 20:
            continue
        print(f"\n{label} :")
        print(f"  γ à R=100   : {gammas[95]:.4f}")
        print(f"  γ à R=400   : {gammas[395]:.4f}")
        print(f"  γ à R=800   : {gammas[795]:.4f}")
        print(f"  cible       : {target:.4f}")
        # tendance : différence entre moyenne début et moyenne fin
        early = np.nanmean(gammas[100:300])
        late = np.nanmean(gammas[700:1000])
        print(f"  <γ> [100,300]  = {early:.4f}")
        print(f"  <γ> [700,1000] = {late:.4f}")
        print(f"  Δ = {late - early:+.4f}   (négatif = γ augmente)")

    print("""
INTERPRÉTATION
--------------
1. Si γ_eff(R) augmente lentement avec R et tend vers la valeur
   asymptotique (3/2 pour l'aire, 4/3 pour le périmètre), alors la
   loi de puissance n'est qu'une approximation valable à grande
   échelle. Les corrections sous-dominantes dominent aux échelles
   intermédiaires.

2. Si γ_eff(R) oscille autour d'une valeur stable, cela signifie
   que la loi de puissance est bien respectée à chaque échelle,
   mais avec un exposant effectif qui dépend du choix de la fenêtre.

3. Si γ_eff(R) décroît vers zéro, alors la déviation ne suit pas
   une loi de puissance du tout : elle est dominée par un biais
   constant (cas des métriques Manhattan/Chebyshev).

4. Le fait que γ_eff dépende de la taille de la fenêtre (L ou c)
   est un diagnostic de multifractalité : il n'y a pas un seul
   exposant, mais un spectre d'exposants.
""")

if __name__ == "__main__":
    main()