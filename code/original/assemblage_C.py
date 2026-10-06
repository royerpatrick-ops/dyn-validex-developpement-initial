"""Assemblage final de C∞ = (6/π²) ∫ a G(a) da :
 - a < 6 : intégration directe (Gauss-Legendre en ln a, tirages uniformes) + queue analytique a → 0 ;
 - a ≥ 6 : G(a) = L/a³ mesuré par échantillonnage ciblé, contribution (6/π²) L / 6."""
import math, numpy as np
from formule_integrale import echantillon
from queue_cible import G_cible

rng = np.random.default_rng(11); k6 = 6 / math.pi ** 2
# 1. le canal ciblé est-il toute la contribution à a = 6 et 8 ?
for a in [6.0, 8.0]:
    c, ar = echantillon(a, 4_000_000, rng)
    gt, et = G_cible(a, 1_000_000, rng, 5, 5)
    print(f"a = {a}: G uniforme = {c.mean():.3e} ± {c.std()/math.sqrt(len(c)):.1e}   G ciblé = {gt:.3e} ± {et:.1e}")
# 2. L = a³ G(a) pour a ≥ 6
Ls = []
for a in [6, 8, 11, 16, 23, 32, 45, 64]:
    g, e = G_cible(float(a), 1_000_000, rng, 5, 5); Ls.append((a ** 3 * g, a ** 3 * e))
L = np.array(Ls); Lm = np.sum(L[:, 0] / L[:, 1] ** 2) / np.sum(1 / L[:, 1] ** 2); Le = 1 / math.sqrt(np.sum(1 / L[:, 1] ** 2))
print("a³G :", np.round(L[:, 0], 4).tolist(), f" -> L = {Lm:.4f} ± {Le:.4f}")
# 3. a < 6 : intégration directe
bornes = np.linspace(-12, math.log(6), 30); xg, wg = np.polynomial.legendre.leggauss(6)
I = 0.0; V = 0.0; IN = 0.0
for u0, u1 in zip(bornes[:-1], bornes[1:]):
    for x, w in zip(xg, wg):
        u = (u0 + u1) / 2 + (u1 - u0) / 2 * x; a = math.exp(u); jac = (u1 - u0) / 2 * w * a * a
        c, ar = echantillon(a, 2_000_000, rng)
        I += jac * c.mean(); V += (jac * c.std() / math.sqrt(len(c))) ** 2; IN += jac * ar.mean()
a0 = math.exp(-12); petit = 8 * math.sqrt(2) / 15 * math.sqrt(a0)
grand = Lm / 6
C = k6 * (I + petit + grand); eC = k6 * math.sqrt(V + (Le / 6) ** 2)
print(f"\nparts : a < e^-12 : {k6*petit:.5f} ; e^-12 < a < 6 : {k6*I:.5f} ± {k6*math.sqrt(V):.5f} ; a > 6 : {k6*grand:.5f} ± {k6*Le/6:.5f}")
print(f"C∞ = {C:.4f} ± {eC:.4f}      (densité de sommets sur a < 6 : {k6*IN:.5f})")
