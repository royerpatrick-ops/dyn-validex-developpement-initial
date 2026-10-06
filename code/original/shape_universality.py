"""
Étape 1C : Test de forme — l'exposant γ est-il universel ?

Génère différentes formes convexes sur un réseau euclidien,
calcule le déficit entre le périmètre continu et le périmètre
de l'enveloppe convexe des points entiers, et ajuste γ.

Dépendances : numpy, matplotlib, scipy
"""

try:
    import numpy as np  # type: ignore[reportMissingImports]
except ImportError as exc:
    raise ImportError(
        "NumPy is required. Install it in the active Python environment with: "
        "python -m pip install numpy"
    ) from exc
try:
    import matplotlib.pyplot as plt  # type: ignore[reportMissingImports]
except ImportError as exc:
    raise ImportError(
        "Matplotlib is required. Install it in the active Python environment with: "
        "python -m pip install matplotlib"
    ) from exc
from scipy.spatial import ConvexHull  # type: ignore[reportMissingImports]
from scipy.ndimage import binary_erosion  # type: ignore[reportMissingImports]
try:
    from scipy.special import ellipe  # type: ignore[reportMissingImports]
except ImportError as exc:
    raise ImportError(
        "SciPy is required. Install it in the active Python environment with: "
        "python -m pip install scipy"
    ) from exc
import time

# ============================================================
# 1. GÉNÉRATION DES FORMES DISCRÈTES
# ============================================================

def grid(R):
    L = int(np.ceil(R)) + 2
    x = np.arange(-L, L + 1)
    return np.meshgrid(x, x)

def gen_circle(R):
    X, Y = grid(R)
    return X*X + Y*Y <= R*R

def gen_ellipse(R, ratio=2.0):
    X, Y = grid(R)
    a, b = R, R / ratio
    return (X/a)**2 + (Y/b)**2 <= 1

def gen_rotated_square(R, angle_deg):
    X, Y = grid(R)
    th = np.deg2rad(angle_deg)
    c, s = np.cos(th), np.sin(th)
    Xr = c*X + s*Y
    Yr = -s*X + c*Y
    # circumradius = R  →  inradius = R/√2
    side = R / np.sqrt(2)
    return (np.abs(Xr) <= side) & (np.abs(Yr) <= side)

def gen_regular_polygon(R, n, angle_deg=0.0):
    X, Y = grid(R)
    th0 = np.deg2rad(angle_deg)
    inside = np.ones(X.shape, dtype=bool)
    r_in = R * np.cos(np.pi/n)   # inradius pour circumradius R
    for k in range(n):
        th = th0 + 2*np.pi*k/n
        nx, ny = np.cos(th), np.sin(th)
        inside &= (nx*X + ny*Y) <= r_in
    return inside

# ============================================================
# 2. PÉRIMÈTRE DE L'ENVELOPPE CONVEXE
# ============================================================

def hull_perimeter(mask):
    """Périmètre de l'enveloppe convexe du masque binaire."""
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
    return float(np.sqrt((d**2).sum(axis=1)).sum())

# ============================================================
# 3. PÉRIMÈTRES CONTINUS (FORMES FERMÉES)
# ============================================================

def P_circle(R):
    return 2 * np.pi * R

def P_ellipse(R, ratio=2.0):
    a, b = R, R / ratio
    e2 = 1 - (b/a)**2
    return 4 * a * ellipe(e2)

def P_rotated_square(R):
    # carré de circumradius R  →  côté = R√2
    return 4 * R * np.sqrt(2)

def P_regular_polygon(R, n):
    return 2 * n * R * np.sin(np.pi/n)

# ============================================================
# 4. PROGRAMME PRINCIPAL
# ============================================================

def main():
    R_max = 500
    R_values = np.arange(10, R_max + 1)
    n_R = len(R_values)

    shapes = [
        ("Cercle",          gen_circle,          {},
         P_circle,          {},              'b'),
        ("Ellipse a/b=2",   gen_ellipse,         {'ratio': 2.0},
         P_ellipse,         {'ratio': 2.0},  'r'),
        ("Ellipse a/b=3",   gen_ellipse,         {'ratio': 3.0},
         P_ellipse,         {'ratio': 3.0},  'm'),
        ("Carré rot. 15°",  gen_rotated_square,  {'angle_deg': 15.0},
         P_rotated_square,  {},              'g'),
        ("Carré rot. 30°",  gen_rotated_square,  {'angle_deg': 30.0},
         P_rotated_square,  {},              'c'),
        ("Pentagone n=5",   gen_regular_polygon, {'n': 5, 'angle_deg': 0.0},
         P_regular_polygon, {'n': 5},        'orange'),
        ("Hexagone n=6",    gen_regular_polygon, {'n': 6, 'angle_deg': 0.0},
         P_regular_polygon, {'n': 6},        'k'),
        ("Polygone n=20",   gen_regular_polygon, {'n': 20, 'angle_deg': 0.0},
         P_regular_polygon, {'n': 20},       'purple'),
        ("Polygone n=50",   gen_regular_polygon, {'n': 50, 'angle_deg': 0.0},
         P_regular_polygon, {'n': 50},       'brown'),
    ]

    print("=" * 70)
    print("TEST DE FORME : UNIVERSALITÉ DE L'EXPOSANT γ")
    print("=" * 70)

    results = []
    t0 = time.time()
    for name, gen, gkw, Pc, pkw, col in shapes:
        print(f"\nForme : {name}")
        Ph_arr = np.empty(n_R)
        Pc_arr = np.empty(n_R)
        for i, R in enumerate(R_values):
            mask = gen(R, **gkw)
            Ph_arr[i] = hull_perimeter(mask)
            Pc_arr[i] = Pc(R, **pkw)
        deficit = np.abs(1.0 - Ph_arr / Pc_arr)
        results.append((name, R_values.astype(float), deficit, col))
        print(f"  Terminé en {time.time()-t0:.1f}s   "
              f"déficit final = {deficit[-1]:.3e}")

    # ============================================================
    # AJUSTEMENTS
    # ============================================================
    print("\n" + "=" * 70)
    print("AJUSTEMENT |1 - P_hull/P_cont| ~ α · R^(-γ)   sur [100, 500]")
    print("=" * 70)

    gammas = {}
    for name, R_arr, deficit, col in results:
        mask = (R_arr >= 100) & (deficit > 1e-15)
        if np.sum(mask) < 10:
            print(f"  {name:<20s} : données insuffisantes")
            gammas[name] = (np.nan, np.nan)
            continue
        logR = np.log(R_arr[mask])
        logD = np.log(deficit[mask])
        A = np.vstack([logR, np.ones_like(logR)]).T
        coef, _, _, _ = np.linalg.lstsq(A, logD, rcond=None)
        gamma = -coef[0]
        alpha = np.exp(coef[1])
        gammas[name] = (alpha, gamma)
        print(f"  {name:<20s} : α = {alpha:.4e}   γ = {gamma:.4f}")

    # ============================================================
    # TRACÉS
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax = axes[0]
    for name, R_arr, deficit, col in results:
        ax.loglog(R_arr, deficit, '.', color=col, markersize=2, label=name)
    ax.set_xlabel('R')
    ax.set_ylabel('|1 − P_hull / P_continu|')
    ax.set_title('Déficit relatif par forme')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[1]
    for name, R_arr, deficit, col in results:
        a, g = gammas[name]
        lab = f"{name} (γ={g:.3f})" if not np.isnan(g) else name
        ax.loglog(R_arr, deficit, '.', color=col, markersize=2, label=lab)
        if not np.isnan(g):
            m = (R_arr >= 100) & (deficit > 1e-15)
            Rf = R_arr[m]
            ax.loglog(Rf, a * Rf**(-g), '-', color=col, linewidth=1.2)
    ax.set_xlabel('R')
    ax.set_ylabel('Déficit relatif')
    ax.set_title('Ajustements γ (pente)')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('shape_universality.png', dpi=150)
    plt.show()

    print("""
INTERPRÉTATION
--------------
1. Si toutes les formes donnent γ ≈ 4/3 : l'exposant est UNIVERSEL
   — propriété purement arithmétique du réseau, indépendante de la
   forme.

2. Si γ ≈ 4/3 pour les formes lisses et γ ≈ 1 pour les polygones
   tournés : DEUX RÉGIMES. Le régime "lisse" est dominé par la
   courbure (Jarník), le régime "polygonal" par les coins
   (corner-cutting O(1)).

3. Si γ varie continûment avec la "rugosité" de la forme (n grand
   pour les polygones réguliers se rapproche du cercle), alors il
   existe un SPECTRE d'exposants qui reflète la structure
   multi-échelle de la géométrie discrète.

4. Conséquence pour la conjecture : l'exposant γ n'est pas une
   constante universelle, mais une FONCTIONNELLE de la forme du
   domaine et de l'observable. Il faut donc écrire :

       π_eff(R, K, O) = π [1 + Σ_k α_k(K, O) · (ℓ_P/R)^{γ_k(K, O)}]

   où K est la forme du domaine et O l'observable.
""")

if __name__ == "__main__":
    main()