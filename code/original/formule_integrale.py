"""
C∞ comme intégrale explicite (somme sur les vecteurs primitifs, à la manière de Huxley).

Modèle : parabole y = -x²/2, réseau affine unimodulaire aléatoire (Haar + translation).
Pour chaque vecteur primitif w, le côté de l'enveloppe de normale ⊥ w se trouve sur la
plus haute droite du réseau parallèle à w qui contient un point sous la parabole.
Par cisaillement (qui préserve la parabole, les écarts verticaux et dx), on se ramène à
w horizontal, de longueur a. Les droites sont alors horizontales, espacées de h = 1/a :
  droite k à la profondeur d_k = δ + k/a,   points en x = s + kβ + j a  (j ∈ ℤ),
avec δ ~ U[0, 1/a), s ~ U[0, a) (translation) et β ~ U[0, a) (forme du réseau, Siegel).
La première droite k dont la corde |x| ≤ √(2 d_k) contient un point porte le côté :
arête si ≥ 2 points (x1 < x2), sinon simple sommet.
Contribution de l'arête au déficit :  c = ∫_{x1}^{x2} (g - g'²/2) dx = d (x2 - x1) - (x2³ - x1³)/3.

Formule de Siegel pour les paires (réseau, vecteur primitif marqué) :
  E[Σ_w F] = (1/ζ(2)) ∫_{R²} dw  E_β[F],  et le changement de variable de cisaillement
  donne, par unité de longueur le long de la parabole :
      C∞ = (6/π²) ∫_0^∞ a · E_{δ,s,β}[ c(a, δ, s, β) ] da
      N  = (6/π²) ∫_0^∞ a · P_{δ,s,β}[ arête ] da      (densité de sommets = densité d'arêtes)
N doit valoir c0/2π ≈ 0,5497 (Balog–Deshouillers).
"""
import math, sys, time
import numpy as np  # pyright: ignore[reportMissingImports]


def echantillon(a, n, rng, kmax=200000):
    """Pour un a donné : n tirages (δ, s, β) -> contribution c et indicatrice d'arête."""
    delta = rng.uniform(0, 1 / a, n); s = rng.uniform(0, a, n); beta = rng.uniform(0, a, n)
    c = np.zeros(n); arete = np.zeros(n, bool)
    actif = np.ones(n, bool); k = 0
    while actif.any() and k < kmax:
        idx = np.where(actif)[0]
        d = delta[idx] + k / a
        r = np.sqrt(2 * d)
        off = np.mod(s[idx] + k * beta[idx], a)          # position du point le plus à gauche ≥ 0
        # points x = off + j a ; ceux dans [-r, r]
        jmin = np.ceil((-r - off) / a); jmax = np.floor((r - off) / a)
        npts = jmax - jmin + 1
        trouve = npts >= 1
        deux = npts >= 2
        sel = idx[trouve & deux]
        if sel.size:
            x1 = off[trouve & deux] + jmin[trouve & deux] * a
            x2 = off[trouve & deux] + jmax[trouve & deux] * a
            dd = d[trouve & deux]
            c[sel] = dd * (x2 - x1) - (x2 ** 3 - x1 ** 3) / 3
            arete[sel] = True
        actif[idx[trouve]] = False
        k += 1
    if actif.any():
        raise RuntimeError(f"a = {a}: {actif.sum()} tirages sans droite trouvée")
    return c, arete


def main(n_par_noeud=200000, graine=3):
    rng = np.random.default_rng(graine)
    # intégration en u = ln a (Gauss-Legendre par morceaux), a de e^-12 à e^2.5
    bornes = np.linspace(-12, math.log(6), 30)   # a ≤ 6 ; au-delà, voir plus bas
    xg, wg = np.polynomial.legendre.leggauss(6)
    I_C = 0.0; I_N = 0.0; V_C = 0.0; V_N = 0.0
    t0 = time.time(); lignes = []
    for u0, u1 in zip(bornes[:-1], bornes[1:]):
        pC = pN = 0.0
        for x, w in zip(xg, wg):
            u = (u0 + u1) / 2 + (u1 - u0) / 2 * x
            a = math.exp(u); jac = (u1 - u0) / 2 * w * a * a      # da = a du, poids a
            c, ar = echantillon(a, n_par_noeud, rng)
            I_C += jac * c.mean(); V_C += (jac * c.std() / math.sqrt(n_par_noeud)) ** 2
            I_N += jac * ar.mean(); V_N += (jac * ar.std() / math.sqrt(n_par_noeud)) ** 2
            pC += jac * c.mean(); pN += jac * ar.mean()
        lignes.append((math.exp(u0), math.exp(u1), pC * 6 / math.pi ** 2, pN * 6 / math.pi ** 2))
        print(f"  a ∈ [{math.exp(u0):9.2e}, {math.exp(u1):6.2f}]  part de C∞ = {pC*6/math.pi**2:.5f}"
              f"   part de N = {pN*6/math.pi**2:.5f}   ({time.time()-t0:.0f} s)", flush=True)
    k = 6 / math.pi ** 2
    a0 = math.exp(bornes[0])
    queue = 8 * math.sqrt(2) / 15 * math.sqrt(a0)   # ∫_0^a0 a·G(a) da avec G ≈ (4√2/15) a^(-3/2)
    print(f"queue analytique a < {a0:.1e} : {k*queue:.5f}")
    I_C += queue
    # a > 6 : l'échantillonnage uniforme de β ne voit presque jamais les rares configurations
    # qui contribuent (deux colonnes de points à ±a/2) ; on les échantillonne directement.
    from queue_cible import G_cible
    Ls = []
    for a in [6.0, 11.0, 23.0, 45.0]:
        g, e = G_cible(a, 200000, rng, 5, 5); Ls.append((a ** 3 * g, a ** 3 * e))
    L = np.array(Ls); Lm = np.sum(L[:, 0] / L[:, 1] ** 2) / np.sum(1 / L[:, 1] ** 2)
    Le = 1 / math.sqrt(np.sum(1 / L[:, 1] ** 2))
    print(f"a > 6 : G(a) ≈ L/a³ avec L = {Lm:.4f} ± {Le:.4f}, part de C∞ = {k*Lm/6:.5f}")
    I_C += Lm / 6; V_C += (Le / 6) ** 2
    print(f"\nN (a ≤ 6, il manque ~0,0001 au-delà) = {k*I_N:.5f} ± {k*math.sqrt(V_N):.5f}   (c0/2π = {3.4536898915/(2*math.pi):.5f})")
    print(f"C∞ = {k*I_C:.5f} ± {k*math.sqrt(V_C):.5f}")
    np.save("formule_integrale.npy", np.array(lignes))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200000)
