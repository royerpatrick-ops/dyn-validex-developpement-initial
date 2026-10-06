"""
Analyse du balayage grands R : départage des deux lois de correction
    (A) C_eff = C_inf + b · R^(-1/3)
    (B) C_eff = C_inf + b / ln R

Données : balayage_translations.csv (R de 1 000 à 30 000, 64 tirages)
        + balayage_grands_R.csv     (R de 30 000 à 10^6, 128 tirages).

Trois tests, pour chaque forme puis pour toutes les formes ensemble :
  1. Prédiction : chaque loi est ajustée sur R ≤ 30 000 seulement, puis
     comparée aux mesures à R > 30 000 (écart en nombre d'erreurs-types).
  2. Ajustement global sur R ≥ 1 000 : χ² réduit de chaque loi.
  3. Pente locale dC/d(ln R) sur R ≥ 100 000, comparée à la pente prédite
     par chaque loi (elle décroît comme R^(-1/3) en (A), comme 1/ln²R en (B)).
"""
import numpy as np  # pyright: ignore[reportMissingImports]
import pandas as pd  # pyright: ignore[reportMissingImports]
from balayage_translations import invariants

LOIS = {"A : R^-1/3": lambda R: R ** (-1 / 3), "B : 1/ln R": lambda R: 1 / np.log(R)}


def charger():
    d1 = pd.read_csv("balayage_translations.csv")
    d1 = d1[d1.R >= 1000].copy()
    # le disque est la superellipse p = 2 dans le premier balayage, l'ellipse 1:1 dans le second
    disque = (d1.forme == "super") & (d1.param == 2.0)
    d1.loc[disque, "forme"] = "ellipse"; d1.loc[disque, "param"] = 1.0
    d2 = pd.read_csv("balayage_grands_R.csv")
    d = pd.concat([d1, d2], ignore_index=True)
    d = d[d.set_index(["forme", "param"]).index.isin(d2.set_index(["forme", "param"]).index.unique())]
    I = {k: invariants(*k)[2] for k in d.groupby(["forme", "param"]).groups}
    d["C"] = d.deficit_abs * d.R ** (1 / 3) / [I[(f, p)] for f, p in zip(d.forme, d.param)]
    g = d.groupby(["forme", "param", "R"]).C.agg(["mean", "sem", "count"]).reset_index()
    return g


def ajuste(R, C, se, X):
    A = np.c_[np.ones_like(R), X(R)]; W = 1 / se
    c, *_ = np.linalg.lstsq(A * W[:, None], C * W, rcond=None)
    chi2 = float(np.sum(((A @ c - C) * W) ** 2))
    cov = np.linalg.inv((A * W[:, None]).T @ (A * W[:, None]))
    return c, chi2, cov


def main():
    g = charger()
    resume = {loi: {"z": [], "chi2": 0.0, "ddl": 0} for loi in LOIS}
    for (f, p), s in g.groupby(["forme", "param"]):
        R, C, se = s.R.values.astype(float), s["mean"].values, s["sem"].values
        print(f"\n=== {f} {p:g} ===")
        print("   R        C_eff     erreur   (tirages)")
        for r, c, e, n in zip(R, C, se, s["count"]):
            if r >= 30000:
                print(f"  {int(r):>8d}  {c:.4f}  ±{e:.4f}   ({n})")
        bas, haut = R <= 30000, R > 30000
        for loi, X in LOIS.items():
            # 1. prédiction depuis R ≤ 30 000
            c, _, cov = ajuste(R[bas], C[bas], se[bas], X)
            A = np.c_[np.ones(haut.sum()), X(R[haut])]
            pred = A @ c; pred_se = np.sqrt(np.einsum("ij,jk,ik->i", A, cov, A))
            z = (C[haut] - pred) / np.sqrt(se[haut] ** 2 + pred_se ** 2)
            resume[loi]["z"].extend(z)
            # 2. ajustement global
            c2, chi2, cov2 = ajuste(R, C, se, X)
            resume[loi]["chi2"] += chi2; resume[loi]["ddl"] += len(R) - 2
            print(f"  {loi}: prédiction à R=10^6 depuis R≤30 000 = {pred[-1]:.4f} ± {pred_se[-1]:.4f}"
                  f" (mesuré {C[-1]:.4f}, écart {z[-1]:+.1f} σ) ;"
                  f" ajustement global C_inf = {c2[0]:.4f} ± {np.sqrt(cov2[0,0]):.4f}, χ²/ddl = {chi2/(len(R)-2):.2f}")
        # 3. pente locale sur R ≥ 100 000
        m = R >= 100000
        if m.sum() >= 3:
            A = np.c_[np.ones(m.sum()), np.log(R[m])]; W = 1 / se[m]
            c, *_ = np.linalg.lstsq(A * W[:, None], C[m] * W, rcond=None)
            cov = np.linalg.inv((A * W[:, None]).T @ (A * W[:, None]))
            Rm = np.exp(np.mean(np.log(R[m])))
            preds = []
            for loi, X in LOIS.items():
                cc, _, _ = ajuste(R, C, se, X)
                h = 1e-4
                preds.append(f"{loi.split(':')[0].strip()} {(cc[1]*(X(Rm*np.exp(h))-X(Rm*np.exp(-h)))/(2*h)):.4f}")
            print(f"  pente dC/d(ln R) mesurée sur R ≥ 10^5 : {c[1]:.4f} ± {np.sqrt(cov[1,1]):.4f}"
                  f"   (prédite : {', '.join(preds)})")
    print("\n=== Bilan, toutes formes ===")
    for loi, r in resume.items():
        z = np.array(r["z"])
        print(f"  {loi}: prédictions à R > 30 000 : écart moyen {z.mean():+.2f} σ, Σz²/n = {np.mean(z**2):.2f} ;"
              f" ajustement global χ²/ddl = {r['chi2']/r['ddl']:.2f}")
    print("""
Lecture : la bonne loi doit donner des écarts de prédiction centrés sur 0
avec Σz²/n ≈ 1, et un χ²/ddl global ≈ 1. Un écart moyen nettement positif
signifie que C_eff monte plus vite que prévu (en faveur de B) ; nettement
négatif, qu'il plafonne plus tôt (en faveur de A).""")


if __name__ == "__main__":
    main()
