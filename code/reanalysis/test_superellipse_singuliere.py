"""
Superellipse |x|^p + |y|^p ≤ 1 avec p < 2 : test de la loi à deux termes aux
grands R, en présence de points de courbure infinie (3 octobre 2026).

Pour 1 < p < 2, la courbure diverge aux quatre points des axes comme |x|^(p−2).
Les intégrales I_4/3 et I_3/2 restent finies si p > 4/3, mais le modèle local
cesse de s'appliquer dans un voisinage de ces points, là où l'échelle locale
κ_a^(−1/3) devient comparable à la distance au point singulier. En coordonnées
réduites, ce voisinage a une largeur ε ~ R^(−2/(p+1)). La part de ∫κ^(4/3) ds
qu'il porte est d'ordre ε^((4p−5)/3), d'où une correction attendue

    C_eff − C_théo ≈ A · R^(−γ),   γ(p) = 2(4p − 5) / (3(p + 1)),

soit γ = 4/15 ≈ 0,267 pour p = 1,5 (une prédiction d'échelle : ni A ni son
signe ne sont prédits). Ce script mesure C_eff de R = 10³ à 10⁶ et ajuste
l'écart à la loi à deux termes.

Utilisation
  python test_superellipse_singuliere.py                 -> p = 1,5, plan complet
  python test_superellipse_singuliere.py --rapide        -> 128 tirages par rayon
  python test_superellipse_singuliere.py --p 1.75        -> autre exposant
  python test_superellipse_singuliere.py --coeurs 4      -> nombre de processus
Reprise automatique : les tirages déjà calculés sont relus dans le fichier JSON.

Dépendances : numpy, scipy, matplotlib ; enveloppe_commune.py dans le même dossier.
"""
import sys, os, json, math, time, argparse
from multiprocessing import Pool
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import numpy as np
from scipy.optimize import curve_fit
import enveloppe_commune as EC

PLAN = [(1000, 1024), (3000, 1024), (10000, 1024), (30000, 1024),
        (100000, 1024), (300000, 512), (1000000, 256)]


def _tirage(args):
    p, R, k = args
    f = EC.superellipse(p)
    rng = np.random.default_rng([2026, int(round(1000 * p)), int(R), k])
    th = rng.uniform(0, np.pi / 2); u, v = rng.uniform(0, 1, 2)
    return k, EC.perimetre_tirage(R, f, th, u, v, n_or=48)


def gamma_pred(p):
    return 2 * (4 * p - 5) / (3 * (p + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=float, default=1.5)
    ap.add_argument("--coeurs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--rapide", action="store_true")
    a = ap.parse_args()
    p = a.p
    plan = [(R, 128 if a.rapide else n) for R, n in PLAN]

    EC.controle()
    P1, I43, I32 = EC.invariants_superellipse(p)
    ratio = I32 / I43
    print(f"Superellipse p = {p} : P = {P1:.8f}, I_4/3 = {I43:.6f}, I_3/2 = {I32:.6f}, ratio = {ratio:.6f}")
    print(f"Exposant prédit de la correction singulière : γ = {gamma_pred(p):.4f}")
    print(f"{a.coeurs} processus\n")

    fichier = f"superellipse_p{p}_singuliere.json"
    donnees = json.load(open(fichier)) if os.path.exists(fichier) else {}

    with Pool(a.coeurs) as pool:
        for R, n in plan:
            cle = str(R)
            faits = {int(k): v for k, v in donnees.get(cle, {}).items()}
            reste = [(p, R, k) for k in range(n) if k not in faits]
            t0 = time.time()
            for i, (k, Ph) in enumerate(pool.imap_unordered(_tirage, reste, chunksize=4)):
                faits[k] = Ph
                if (i + 1) % 128 == 0 or i + 1 == len(reste):
                    donnees[cle] = {str(kk): vv for kk, vv in faits.items()}
                    json.dump(donnees, open(fichier, "w"))
                    print(f"  R = {R:>7} : {len(faits)}/{n} tirages ({time.time() - t0:.0f} s)", flush=True)

    # ------------------------------------------------------------------ bilan
    rows = []
    print(f"\n{'R':>8} {'tirages':>8} {'C_eff':>17} {'C_théo':>8} {'ΔC':>9} {'résidu':>8}")
    for R, n in plan:
        Ph = np.array([donnees[str(R)][str(k)] for k in range(n)])
        fac = R ** (1 / 3) / I43
        ce = (P1 * R - Ph.mean()) * fac
        se = Ph.std(ddof=1) / math.sqrt(n) * fac
        ct = EC.C_theo(R, ratio)
        rows.append((R, ce, se, ct))
        print(f"{R:>8} {n:>8} {ce:8.4f} ± {se:.4f} {ct:8.4f} {ce - ct:+9.4f} {(ce - ct) / se:+7.2f}σ")

    R_ = np.array([r[0] for r in rows], float)
    d_ = np.array([r[1] - r[3] for r in rows])
    s_ = np.array([r[2] for r in rows])
    g0 = gamma_pred(p)
    chi_nul = float(np.sum((d_ / s_) ** 2))
    w = 1 / s_ ** 2
    x0 = R_ ** (-g0)
    A0 = float(np.sum(w * d_ * x0) / np.sum(w * x0 * x0))
    eA0 = float(1 / math.sqrt(np.sum(w * x0 * x0)))
    chi_0 = float(np.sum(((d_ - A0 * x0) / s_) ** 2))
    try:
        (A1, g1), cov = curve_fit(lambda R, A, g: A * R ** (-g), R_, d_, p0=[A0, g0], sigma=s_,
                                  absolute_sigma=True, maxfev=20000)
        eA1, eg1 = np.sqrt(np.diag(cov))
        chi_1 = float(np.sum(((d_ - A1 * R_ ** (-g1)) / s_) ** 2))
    except Exception:
        A1 = g1 = eA1 = eg1 = chi_1 = float("nan")
    n = len(rows)
    print(f"\nÉcart nul (loi à deux termes seule)     : χ² = {chi_nul:.1f} pour {n} points")
    print(f"ΔC = A·R^(−{g0:.3f}) (exposant prédit)    : A = {A0:+.4f} ± {eA0:.4f}, χ² = {chi_0:.1f} ({n - 1} ddl)")
    print(f"ΔC = A·R^(−γ), γ libre                  : A = {A1:+.4f} ± {eA1:.4f}, γ = {g1:.3f} ± {eg1:.3f}, "
          f"χ² = {chi_1:.1f} ({n - 2} ddl)")

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    Rt = np.geomspace(800, 1.3e6, 200)
    ax[0].errorbar(R_, [r[1] for r in rows], yerr=s_, fmt='o', color='C3', capsize=4, label=f'p = {p} (mesuré)')
    ax[0].plot(Rt, EC.C_INF - EC.B_CONST * ratio * Rt ** (-1 / 6), '--', color='k', alpha=0.6, label='loi à deux termes')
    ax[0].plot(Rt, EC.C_INF - EC.B_CONST * ratio * Rt ** (-1 / 6) + A0 * Rt ** (-g0), '-', color='C0', alpha=0.7,
               label=f'+ A·R^(−{g0:.3f})')
    ax[0].axhline(EC.C_INF, color='k', linestyle=':', alpha=0.4)
    ax[0].set_xscale('log'); ax[0].set_xlabel('R'); ax[0].set_ylabel('C_eff(R)')
    ax[0].set_title(f'Superellipse p = {p}'); ax[0].legend(); ax[0].grid(True, alpha=0.3, which='both')
    ax[1].errorbar(R_, d_, yerr=s_, fmt='o', color='C3', capsize=4, label='C_eff − C_théo')
    ax[1].plot(Rt, A0 * Rt ** (-g0), '-', color='C0', label=f'A·R^(−{g0:.3f}) (exposant prédit)')
    if np.isfinite(g1):
        ax[1].plot(Rt, A1 * Rt ** (-g1), '--', color='C2', label=f'γ libre = {g1:.3f}')
    ax[1].axhline(0, color='k', linewidth=0.6)
    ax[1].set_xscale('log'); ax[1].set_xlabel('R'); ax[1].set_ylabel('ΔC')
    ax[1].set_title('Écart à la loi à deux termes'); ax[1].legend(); ax[1].grid(True, alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(f'superellipse_p{p}_singuliere.png', dpi=150)
    print(f"\nFigure : superellipse_p{p}_singuliere.png")
    plt.show()


if __name__ == "__main__":
    main()
