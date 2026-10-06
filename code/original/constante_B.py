"""Prédiction analytique de B (terme en R^(-1/2) ∫κ^(3/2) ds) à partir de la fonction zêta d'Epstein de Z²."""
import math
import numpy as np  # pyright: ignore[reportMissingImports]
import mpmath as mp  # pyright: ignore[reportMissingImports]
mp.mp.dps = 20
z = mp.zeta(0.75); b = mp.dirichlet(0.75, [0, 1, 0, -1])
Z = 4 * z * b; Zp = Z / mp.zeta(1.5)
g0 = 4 * mp.sqrt(2) / 15
I = mp.quad(lambda u: (1 + u * u) ** (-0.75), [0, mp.inf])
K = g0  # une arête de direction w0 n'apparaît qu'au point où la tangente est parallèle à w0
B = -K * Zp / (2 * mp.pi)
print(f"ζ(3/4) = {float(z):.6f}, β(3/4) = {float(b):.6f}, Z(3/2) = 4ζ(3/4)β(3/4) = {float(Z):.6f}")
print(f"Z_prim(3/2) = Z(3/2)/ζ(3/2) = {float(Zp):.6f}")
print(f"K = {float(K):.6f}   ->   B prédit = {float(B):.4f}")
# contrôle numérique de la régularisation (somme sur les vecteurs primitifs moins intégrale)
X = 1500
p = np.arange(-X, X + 1)
tot = 0.0
for q in range(0, X + 1):
    pp = p if q > 0 else p[p > 0]
    r2 = pp * pp + q * q
    m = (r2 <= X * X) & (np.gcd(pp, q) == 1)
    tot += 2 * np.sum(r2[m].astype(float) ** -0.75)
print(f"contrôle : Σ_prim(|w|<{X}) |w|^(-3/2) − (24/π)·X^(1/2) = {tot - 24/math.pi*X**0.5:.4f}   (attendu {float(Zp):.4f})")
