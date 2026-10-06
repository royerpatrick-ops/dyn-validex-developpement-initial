"""Contribution des grands a à B : arêtes à deux colonnes, dont le vecteur partenaire v ∈ Z² (det(w0,v)=1)
fixe m = λ|v|, λ = a²ε/2. Rapport discret/continu de la probabilité d'arête :
   ρ(λ) = (πλ/8) Σ_{v primitif, |v|<1/λ} (1 − λ²|v|²)/|v|      (→ 1 quand λ → 0)
   B_grand = −(6/π²)(1/9)(√2/4) ∫_0^∞ λ^(−3/2) (ρ(λ) − 1) dλ
"""
import math, numpy as np
from scipy.integrate import quad
N = 4000
p = np.arange(-N, N+1); P, Q = np.meshgrid(p, p)
m = (np.gcd(P, Q) == 1)
r = np.sort(np.sqrt((P[m]**2 + Q[m]**2).astype(float)))
r = r[r <= N]
cum_inv = np.concatenate([[0], np.cumsum(1/r)]); cum_r = np.concatenate([[0], np.cumsum(r)])
def rho(lam):
    k = np.searchsorted(r, 1/lam, side="left")
    return math.pi*lam/8 * (cum_inv[k] - lam**2 * cum_r[k])
lam_min = 2/N   # en deçà, ρ ≈ 1 (continuum) à mieux que 1/N près
pts = sorted(set([1/x for x in np.unique(r[r < 1/lam_min])]))
pts = [x for x in pts if lam_min < x < 1.0] + [1.0]
I = 0.0; prev = lam_min
for x in pts:
    I += quad(lambda l: l**-1.5*(rho(l)-1), prev, x, limit=50)[0]; prev = x
I += quad(lambda l: -l**-1.5, 1.0, np.inf)[0]
# correction de la partie λ < lam_min (ρ−1 = O(λ^(1/2))·oscillant, de moyenne nulle) : négligée
Bg = -(6/math.pi**2)*(1/9)*(math.sqrt(2)/4)*I
Bp = 4*math.sqrt(2)/15*3.857623/(2*math.pi)
print(f"∫ λ^(-3/2)(ρ−1) dλ = {I:.4f}   (modèle « trou |v|<1 » uniforme : −4.8)")
print(f"B (petits a, zêta d'Epstein) = {Bp:.4f}")
print(f"B (grands a, partenaires entiers) = {Bg:.4f}")
print(f"B total = {Bp+Bg:.4f}")
for l in [0.05, 0.1, 0.3, 0.5, 0.8, 0.99]: print(f"   ρ({l}) = {rho(l):.3f}")
