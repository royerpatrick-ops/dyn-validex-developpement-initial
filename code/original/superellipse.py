"""
Étape 1E : Balayage continu en p — superellipses.
Version simplifiée.
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
from scipy.spatial import ConvexHull  # type: ignore
from scipy.ndimage import binary_erosion  # type: ignore
import time

# Compatibilité NumPy 2.0
_integrate = getattr(np, "trapezoid", None) or np.trapz


def superellipse_mask(R, p, theta0=0.3):
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


def perimeter_continuous(R, p, N_theta=8000):
    if abs(p - 1.0) < 1e-9:
        return 4 * np.sqrt(2) * R
    if p > 1e8:
        return 8 * R
    q = p / (p - 1)
    theta = np.linspace(0, np.pi, N_theta)
    integrand = (np.abs(np.cos(theta)) ** q
                 + np.abs(np.sin(theta)) ** q) ** (1.0 / q)
    return 2 * R * _integrate(integrand, theta)


def fit_gamma(R_arr, deficit):
    mask = deficit > 1e-14
    if np.sum(mask) < 10:
        return np.nan
    logR = np.log(R_arr[mask])
    logD = np.log(deficit[mask])
    A = np.vstack([logR, np.ones_like(logR)]).T
    coef, _, _, _ = np.linalg.lstsq(A, logD, rcond=None)
    return -coef[0]


def main():
    p_values = np.concatenate([
        np.arange(1.00, 1.50, 0.05),
        np.arange(1.50, 3.00, 0.10),
        np.arange(3.00, 5.00, 0.20),
        np.arange(5.00, 10.01, 0.50),
    ])
    R_values = np.arange(40, 301, dtype=float)

    windows = {
        '40-120':  (R_values >= 40)  & (R_values <= 120),
        '80-200':  (R_values >= 80)  & (R_values <= 200),
        '150-300': (R_values >= 150),
        '80-300':  (R_values >= 80),
    }

    gammas = {name: np.full(len(p_values), np.nan) for name in windows}
    deficits_examples = {}

    print("=" * 70)
    print("BALAYAGE CONTINU EN p -- SUPERELLIPSES")
    print("=" * 70)
    print(f"p de {p_values[0]:.2f} a {p_values[-1]:.2f} "
          f"({len(p_values)} valeurs)")
    print(f"R de {int(R_values[0])} a {int(R_values[-1])}")

    t0 = time.time()
    for i, p in enumerate(p_values):
        defs = np.empty(len(R_values))
        for j, R in enumerate(R_values):
            mask = superellipse_mask(R, p)
            Ph = hull_perimeter(mask)
            Pc = perimeter_continuous(R, p)
            defs[j] = abs(1.0 - Ph / Pc) if Pc > 0 else 0.0
        for name, mk in windows.items():
            gammas[name][i] = fit_gamma(R_values[mk], defs[mk])

        if i % 3 == 0 or p in (1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0):
            deficits_examples[p] = (R_values.copy(), defs.copy())

        if (i + 1) % 3 == 0 or i == 0:
            g = gammas['80-300'][i]
            print(f"  p = {p:5.2f}   gamma[80-300] = {g:.4f}   "
                  f"({time.time() - t0:.1f}s)")

    print("\n" + "=" * 70)
    print("gamma(p) PAR FENETRE")
    print("=" * 70)
    header = f"{'p':>6}"
    for name in windows:
        header += f"  {'g[' + name + ']':>14}"
    print(header)
    print("-" * len(header))
    for i, p in enumerate(p_values):
        if (abs(p - round(p, 1)) < 0.01
                or p in (2.0, 3.0, 5.0, 8.0, 10.0)):
            row = f"{p:>6.2f}"
            for name in windows:
                row += f"  {gammas[name][i]:>14.4f}"
            print(row)

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    ax = axes[0, 0]
    styles = {'40-120': 'o', '80-200': 's', '150-300': '^', '80-300': 'd'}
    for name in windows:
        ax.plot(p_values, gammas[name], styles[name] + '-',
                markersize=4, linewidth=0.9, label=f'fenetre {name}')
    ax.axhline(1.0, color='r', linestyle='--', linewidth=1.2,
               label='gamma = 1 (coins)')
    ax.axhline(4 / 3, color='b', linestyle='--', linewidth=1.2,
               label='gamma = 4/3 (Jarnik)')
    ax.axhline(1.27, color='g', linestyle=':', linewidth=1.2,
               label='gamma ~ 1.27 (cercle mesure)')
    ax.set_xlabel('p')
    ax.set_ylabel('gamma(p)')
    ax.set_title('Transition coins -> courbure (superellipses tournees)')
    ax.set_ylim(0.85, 1.45)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    mask_zoom = p_values <= 2.5
    for name in windows:
        ax.plot(p_values[mask_zoom], gammas[name][mask_zoom],
                styles[name] + '-', markersize=5, linewidth=1.0,
                label=f'fenetre {name}')
    ax.axhline(1.0, color='r', linestyle='--', linewidth=1.2)
    ax.axhline(4 / 3, color='b', linestyle='--', linewidth=1.2)
    ax.set_xlabel('p')
    ax.set_ylabel('gamma(p)')
    ax.set_title('Zoom p dans [1, 2.5]')
    ax.set_ylim(0.9, 1.4)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    for name in windows:
        g = gammas[name]
        valid = ~np.isnan(g)
        if valid.sum() < 5:
            continue
        dg = np.gradient(g[valid], p_values[valid])
        ax.plot(p_values[valid], dg, '-', linewidth=1.0,
                label=f'fenetre {name}')
    ax.axhline(0, color='k', linewidth=0.5)
    ax.set_xlabel('p')
    ax.set_ylabel('dgamma/dp')
    ax.set_title('Derivee dgamma/dp')
    ax.set_xlim(1, 5)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    colors = plt.cm.viridis(np.linspace(0, 1, len(deficits_examples)))
    for k, (p, (R_arr, defs)) in enumerate(sorted(deficits_examples.items())):
        ax.loglog(R_arr, defs, '-', color=colors[k],
                  linewidth=1.0, label=f'p = {p:.2f}')
    ax.set_xlabel('R')
    ax.set_ylabel('Deficit relatif')
    ax.set_title('Deficit pour quelques superellipses')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_ylim(1e-6, 1)

    plt.tight_layout()
    plt.savefig('superellipse_gamma.png', dpi=150)
    plt.show()

    print("\n" + "=" * 70)
    print("ANALYSE DE LA TRANSITION")
    print("=" * 70)
    for name in windows:
        g = gammas[name]
        valid = ~np.isnan(g)
        p_c = p_values[valid]
        g_c = g[valid]
        if len(g_c) < 5:
            print(f"  fenetre {name:8s} : donnees insuffisantes")
            continue
        i_max = np.argmax(g_c)
        i_p2 = np.argmin(np.abs(p_c - 2))
        print(f"  fenetre {name:8s} : "
              f"gamma_max = {g_c[i_max]:.4f} a p = {p_c[i_max]:.2f}   "
              f"gamma(p=2) = {g_c[i_p2]:.4f}")


if __name__ == "__main__":
    main()