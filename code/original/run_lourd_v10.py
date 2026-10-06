"""
Extension à R = 30 000 — version v10 (rapide + robuste).

Méthode : bisection continue pour converger vers y_max/y_min flottants,
puis test des 3 entiers voisins vers l'extérieur pour capturer les
points entiers sur la frontière. O(R) par orientation au lieu de O(R²).

Validation croisée avec grille pleine pour R = 100 (err attendue : 0).

Run lourd : N_orient = 16, 6 formes, R ∈ {300, 1000, 3000, 10000, 30000}.

Dépendances : numpy, scipy, matplotlib
"""

import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import numpy as np
from scipy.spatial import ConvexHull
from scipy.special import ellipe
import matplotlib.pyplot as plt
import time

TOL = 1e-9


def make_shape(kind, param, tol=TOL):
    if kind == "superellipse":
        p = param
        def f(xn, yn):
            return (np.abs(xn) ** p + np.abs(yn) ** p) <= 1.0 + tol
        return f
    elif kind == "ellipse":
        r = param
        def f(xn, yn):
            return (xn ** 2 + (yn * r) ** 2) <= 1.0 + tol
        return f
    raise ValueError(kind)


# ============================================================
# BORNES ENTIÈRES PAR COLONNE (MÉTHODE RAPIDE)
# ============================================================

def column_integer_bounds(xs, R, shape_func, theta0,
                          N_y=400, n_bisect=45):
    """
    Pour chaque x, retourne (y_bot, y_top) entiers.
    Méthode : bisection flottante + test des voisins.
    """
    c, s = np.cos(theta0), np.sin(theta0)
    y_bound = R * np.sqrt(2) + 2

    # --- Passe initiale (linspace + 0 pour garantir un point inside) ---
    ys = np.unique(np.concatenate([
        np.linspace(-y_bound, y_bound, N_y),
        np.array([0.0])
    ]))
    N_y_eff = len(ys)

    X = xs[:, None]
    Y = ys[None, :]
    xp = c * X + s * Y
    yp = -s * X + c * Y
    inside = shape_func(xp / R, yp / R)

    J = np.arange(N_y_eff)
    j_max_raw = np.where(inside, J, -1).max(axis=1)
    j_min_raw = np.where(inside, J, N_y_eff).min(axis=1)
    valid = (j_max_raw >= 0) & (j_min_raw < N_y_eff)

    j_max = np.clip(j_max_raw, 0, N_y_eff - 1)
    j_min = np.clip(j_min_raw, 0, N_y_eff - 1)

    # Encadrement y_max
    j_hi = np.minimum(j_max + 1, N_y_eff - 1)
    y_max_lo = np.where(valid, ys[j_max], -y_bound)
    y_max_hi = np.where(valid, ys[j_hi],  y_bound)
    extend_top = (j_max_raw == N_y_eff - 1) & valid
    y_max_hi = np.where(extend_top, y_bound * 2, y_max_hi)

    # Encadrement y_min
    j_lo = np.maximum(j_min - 1, 0)
    y_min_lo = np.where(valid, ys[j_lo], -y_bound)
    y_min_hi = np.where(valid, ys[j_min],  y_bound)
    extend_bot = (j_min_raw == 0) & valid
    y_min_lo = np.where(extend_bot, -y_bound * 2, y_min_lo)

    # --- Bisection flottante ---
    for _ in range(n_bisect):
        y_mid = 0.5 * (y_max_lo + y_max_hi)
        xp_m = c * xs + s * y_mid
        yp_m = -s * xs + c * y_mid
        in_m = shape_func(xp_m / R, yp_m / R)
        y_max_lo = np.where(in_m & valid, y_mid, y_max_lo)
        y_max_hi = np.where(~in_m & valid, y_mid, y_max_hi)

        y_mid = 0.5 * (y_min_lo + y_min_hi)
        xp_m = c * xs + s * y_mid
        yp_m = -s * xs + c * y_mid
        in_m = shape_func(xp_m / R, yp_m / R)
        y_min_lo = np.where(~in_m & valid, y_mid, y_min_lo)
        y_min_hi = np.where(in_m & valid, y_mid, y_min_hi)

    y_max_cont = 0.5 * (y_max_lo + y_max_hi)
    y_min_cont = 0.5 * (y_min_lo + y_min_hi)

    # --- Test des entiers voisins ---
    y_top = np.floor(y_max_cont)
    for _ in range(3):
        y_test = y_top + 1
        xp_t = c * xs + s * y_test
        yp_t = -s * xs + c * y_test
        in_t = shape_func(xp_t / R, yp_t / R)
        y_top = np.where(in_t & valid, y_test, y_top)

    y_bot = np.ceil(y_min_cont)
    for _ in range(3):
        y_test = y_bot - 1
        xp_b = c * xs + s * y_test
        yp_b = -s * xs + c * y_test
        in_b = shape_func(xp_b / R, yp_b / R)
        y_bot = np.where(in_b & valid, y_test, y_bot)

    y_top = np.where(valid, y_top, np.nan)
    y_bot = np.where(valid, y_bot, np.nan)
    return y_bot, y_top


# ============================================================
# PÉRIMÈTRE DE L'ENVELOPPE
# ============================================================

def hull_perimeter(xs, y_bot, y_top):
    valid = ~np.isnan(y_top) & ~np.isnan(y_bot)
    xs_v = xs[valid]
    if len(xs_v) < 2:
        return 0.0
    pts = np.vstack([
        np.column_stack([xs_v, y_top[valid]]),
        np.column_stack([xs_v, y_bot[valid]])
    ])
    pts = np.unique(pts, axis=0)
    if len(pts) < 3:
        return 0.0
    hull = ConvexHull(pts)
    v = pts[hull.vertices]
    d = np.diff(np.vstack([v, v[:1]]), axis=0)
    return float(np.sqrt((d ** 2).sum(axis=1)).sum())


# ============================================================
# PÉRIMÈTRES CONTINUS ET INVARIANTS
# ============================================================

def P_cont_super(R, p, N=20000):
    if abs(p - 2.0) < 1e-10:
        return 2 * np.pi * R
    q = p / (p - 1)
    t = np.linspace(0, np.pi, N)
    f = (np.abs(np.cos(t)) ** q + np.abs(np.sin(t)) ** q) ** (1.0 / q)
    return 2 * R * np.trapezoid(f, t)


def P_cont_ellipse(R, ratio):
    a, b = R, R / ratio
    return 4 * a * ellipe(1 - (b / a) ** 2)


def invariant_super(p, N=1_000_000):
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    ct, st = np.cos(t), np.sin(t)
    x = np.abs(ct) ** (2.0 / p) * np.sign(ct)
    y = np.abs(st) ** (2.0 / p) * np.sign(st)
    dt = t[1] - t[0]
    xp = np.gradient(x, dt); yp = np.gradient(y, dt)
    xpp = np.gradient(xp, dt); ypp = np.gradient(yp, dt)
    k = np.abs(xp * ypp - yp * xpp) / (xp ** 2 + yp ** 2) ** 1.5
    ds = np.sqrt(xp ** 2 + yp ** 2) * dt
    return float(np.sum(k ** (4.0 / 3.0) * ds))


def invariant_ellipse(ratio, N=1_000_000):
    a, b = 1.0, 1.0 / ratio
    t = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x = a * np.cos(t); y = b * np.sin(t)
    dt = t[1] - t[0]
    xp = -a * np.sin(t); yp = b * np.cos(t)
    xpp = -a * np.cos(t); ypp = -b * np.sin(t)
    k = np.abs(xp * ypp - yp * xpp) / (xp ** 2 + yp ** 2) ** 1.5
    ds = np.sqrt(xp ** 2 + yp ** 2) * dt
    return float(np.sum(k ** (4.0 / 3.0) * ds))


# ============================================================
# VALIDATION
# ============================================================

def validate_full_grid(kind, param, R=100, theta0=0.3):
    shape = make_shape(kind, param)
    R_int = int(R)
    x_bound = int(np.ceil(R_int * np.sqrt(2))) + 2
    y_bound = x_bound

    xs = np.arange(-x_bound, x_bound + 1, dtype=float)
    y_bot, y_top = column_integer_bounds(xs, float(R), shape, theta0)
    P_col = hull_perimeter(xs, y_bot, y_top)

    Xg, Yg = np.meshgrid(np.arange(-x_bound, x_bound + 1),
                         np.arange(-y_bound, y_bound + 1))
    c, s = np.cos(theta0), np.sin(theta0)
    xp = c * Xg + s * Yg
    yp = -s * Xg + c * Yg
    inside = shape(xp / R, yp / R)
    pts = np.column_stack([Xg[inside], Yg[inside]]).astype(float)
    if len(pts) < 3:
        return None
    hull = ConvexHull(pts)
    v = pts[hull.vertices]
    d = np.diff(np.vstack([v, v[:1]]), axis=0)
    P_full = float(np.sqrt((d ** 2).sum(axis=1)).sum())

    return P_col, P_full, abs(P_col - P_full) / P_full


def run_validation():
    print("Validation croisée (colonnes vs grille pleine, R = 100) :")
    all_ok = True
    for kind, param, label in [
        ("superellipse", 2.0, "cercle p=2"),
        ("superellipse", 1.5, "super p=1.5"),
        ("superellipse", 2.5, "super p=2.5"),
        ("superellipse", 5.0, "super p=5"),
        ("superellipse", 10.0, "super p=10"),
        ("ellipse", 3.0, "ellipse a/b=3"),
    ]:
        res = validate_full_grid(kind, param, R=100, theta0=0.3)
        if res:
            Pc, Pf, err = res
            status = "OK" if err < 1e-12 else "ÉCHEC"
            if err >= 1e-12:
                all_ok = False
            print(f"  {label:>14s}  P_col={Pc:.6f}  P_full={Pf:.6f}  "
                  f"err={err:.3e}  [{status}]")
    print()
    if all_ok:
        print(">>> Validation réussie. <<<\n")
    else:
        print(">>> ATTENTION : erreur résiduelle. <<<\n")
    return all_ok


# ============================================================
# PROGRAMME PRINCIPAL
# ============================================================

def main():
    if not run_validation():
        print("Interruption.")
        return

    shapes = [
        ("cercle_p2",    "superellipse", 2.0),
        ("super_p1.5",   "superellipse", 1.5),
        ("super_p2.5",   "superellipse", 2.5),
        ("super_p5",     "superellipse", 5.0),
        ("super_p10",    "superellipse", 10.0),
        ("ellipse_a3",   "ellipse",      3.0),
    ]
    R_values = np.array([300, 1000, 3000, 10000, 30000])
    N_orient = 16
    np.random.seed(12345)

    all_results = {}

    for name, kind, param in shapes:
        print(f"\n=== {name} ===")
        I_K = (invariant_super(param) if kind == "superellipse"
               else invariant_ellipse(param))
        print(f"  I(K) = ∫κ^(4/3) ds = {I_K:.6f}")

        shape = make_shape(kind, param)
        results = []
        t0 = time.time()

        for R in R_values:
            R_int = int(R)
            x_bound = int(np.ceil(R_int * np.sqrt(2))) + 2
            xs = np.arange(-x_bound, x_bound + 1, dtype=float)

            P_cont = (P_cont_super(R_int, param)
                      if kind == "superellipse"
                      else P_cont_ellipse(R_int, param))

            thetas = np.random.uniform(0, np.pi, N_orient)
            Ph_list = []
            for theta in thetas:
                y_bot, y_top = column_integer_bounds(
                    xs, float(R_int), shape, theta)
                Ph_list.append(hull_perimeter(xs, y_bot, y_top))
            P_hull = float(np.mean(Ph_list))

            DeltaP = P_cont - P_hull
            C_eff = DeltaP * (R_int ** (1.0 / 3.0)) / I_K
            results.append((R_int, P_cont, P_hull, DeltaP, C_eff))
            print(f"  R={R_int:>6d}  ΔP={DeltaP:>12.5f}  "
                  f"C_eff={C_eff:.4f}  ({time.time()-t0:.1f}s)")

        all_results[name] = results

    # --- Figure ---
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {
        'cercle_p2':  'b', 'super_p1.5': 'c',
        'super_p2.5': 'orange', 'super_p5': 'g',
        'super_p10':  'purple', 'ellipse_a3': 'r'
    }
    for name, results in all_results.items():
        R_arr = np.array([r[0] for r in results])
        C_arr = np.array([r[4] for r in results])
        ax.semilogx(R_arr, C_arr, 'o-', label=name,
                    color=colors.get(name, None), markersize=8)
    ax.axhline(0.646, color='k', linestyle='--', alpha=0.5,
               label='référence (0.646)')
    ax.set_xlabel('R')
    ax.set_ylabel('C_eff(R)')
    ax.set_title('C_eff — v10, N_orient = 16, 6 formes')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig('C_eff_v10.png', dpi=150)
    plt.show()

    # --- Analyse ---
    print("\n" + "=" * 70)
    print("ANALYSE STATISTIQUE")
    print("=" * 70)
    print(f"\n{'Forme':>14s}  {'C_eff(R=30000)':>15s}  {'pente':>10s}")
    print("-" * 45)
    C_30000 = {}
    for name, results in all_results.items():
        R_arr = np.array([r[0] for r in results])
        C_arr = np.array([r[4] for r in results])
        slope = np.polyfit(np.log(R_arr[-3:]), C_arr[-3:], 1)[0]
        C_30000[name] = C_arr[-1]
        print(f"{name:>14s}  {C_arr[-1]:>15.4f}  {slope:>+10.5f}")

    C_vals = np.array(list(C_30000.values()))
    print(f"\nMoyenne   : {C_vals.mean():.4f}")
    print(f"Écart-type: {C_vals.std():.4f}")
    print(f"Dispersion: {100 * C_vals.std() / C_vals.mean():.2f} %")

    with open("resultats_v10.txt", "w") as f:
        f.write("C_eff v10 — N_orient=16, 6 formes\n")
        f.write("=" * 60 + "\n")
        for name, results in all_results.items():
            f.write(f"\n{name}:\n")
            for R, Pc, Ph, DP, C in results:
                f.write(f"  R={R:>6d}  P_cont={Pc:.4f}  "
                        f"P_hull={Ph:.4f}  DeltaP={DP:.5f}  "
                        f"C_eff={C:.4f}\n")


if __name__ == "__main__":
    main()