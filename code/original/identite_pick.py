"""Test de l'identité de Pick dans le modèle local :
   ⟨g⟩ (déficit d'aire par unité de longueur) = n_∂ / 2,
   n_∂ = densité des points du réseau situés sur le bord de l'enveloppe.
Décomposition : C∞ = ⟨g⟩ − ½⟨g'²⟩."""
import math, sys, numpy as np
from queue_cible import G_cible

def ech(a, n, rng):
    delta = rng.uniform(0, 1/a, n); s = rng.uniform(0, a, n); beta = rng.uniform(0, a, n)
    D = np.zeros(n); Gm = np.zeros(n); Nb = np.zeros(n)
    act = np.arange(n); k = 0
    while act.size:
        d = delta[act] + k/a; r = np.sqrt(2*d); off = np.mod(s[act] + k*beta[act], a)
        jmin = np.ceil((-r-off)/a); jmax = np.floor((r-off)/a); npt = jmax - jmin + 1
        tr = npt >= 1; dx = npt >= 2; sel = act[dx]
        x1 = off[dx] + jmin[dx]*a; x2 = off[dx] + jmax[dx]*a; dd = d[dx]
        D[sel] = dd*(x2-x1) - (x2**3-x1**3)/3
        Gm[sel] = dd*(x2-x1) - (x2**3-x1**3)/6
        Nb[sel] = npt[dx] - 1                      # segments = points du bord comptés une fois par arête
        act = act[~tr]; k += 1
    return D, Gm, Nb

rng = np.random.default_rng(41); k6 = 6/math.pi**2
bornes = np.linspace(-12, math.log(6), 30); xg, wg = np.polynomial.legendre.leggauss(6)
S = np.zeros(3); V = np.zeros(3); n = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
for u0, u1 in zip(bornes[:-1], bornes[1:]):
    for x, w in zip(xg, wg):
        u = (u0+u1)/2 + (u1-u0)/2*x; a = math.exp(u); jac = (u1-u0)/2*w*a*a
        res = ech(a, n, rng)
        for i, r in enumerate(res): S[i] += jac*r.mean(); V[i] += (jac*r.std()/math.sqrt(n))**2
a0 = math.exp(-12)
# queues analytiques (a → 0) : D ≈ (4√2/15) a^-3/2, ∫g ≈ (8√2/15)·... , n ≈ (4√2/3) a^-3/2
tail_small = np.array([4*math.sqrt(2)/15, 0.0, 4*math.sqrt(2)/3*2/3])
# ∫g sur la première droite : d·2r − (2r³)/6 = 2√2 d^{3/2}(1 − 1/3) = (4√2/3) d^{3/2} → moyenne (4√2/3)(2/5) a^-3/2
tail_small[1] = 4*math.sqrt(2)/3*2/5
S += tail_small*2*math.sqrt(a0)
# grands a (> 6) : arêtes à deux colonnes, n = 1, c ≈ a³/24, ∫g ≈ a³/12 ... évalué par échantillonnage ciblé pour c ;
# pour ∫g et n on utilise les rapports asymptotiques exacts : ∫g/c = 2 et n·(a³/24)/c = a³/24 → contribution Σ ∫ a·P da
L = 1/9
S[0] += L/6; S[1] += 2*L/6; S[2] += (8/3)/(4*6**4)   # n : ∫_6^∞ a·P(arête) da, P = 4J/a⁶ = (8/3)/a⁶
print(f"C∞      = {k6*S[0]:.5f} ± {k6*math.sqrt(V[0]):.5f}")
print(f"⟨g⟩     = {k6*S[1]:.5f} ± {k6*math.sqrt(V[1]):.5f}")
print(f"n_∂/2   = {k6*S[2]/2:.5f} ± {k6*math.sqrt(V[2])/2:.5f}    (Pick : doit égaler ⟨g⟩)")
print(f"½⟨g'²⟩  = {k6*(S[1]-S[0]):.5f}")
