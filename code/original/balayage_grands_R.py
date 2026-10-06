"""
Balayage grands R : C_eff pour R de 30 000 à 1 000 000, afin de départager
deux lois de correction de la dérive de C_eff :
    (A) C_eff(R) = C_inf + b · R^(-1/3)     -> C_inf ≈ 0,68
    (B) C_eff(R) = C_inf + b / ln R         -> C_inf ≈ 0,75

Réutilise les fonctions validées de balayage_translations.py (même méthode :
colonnes exactes sur le réseau décalé Z² - (u, v), rotation + translation
aléatoires, graines déterministes).

Formes : disque, ellipses 2:1 et 5:1 (forme fermée, rapides),
superellipse p = 3 (plus lente : ~20 s par tirage à R = 10^6 sur un cœur).

Usage :
    python balayage_grands_R.py                 # tout, tous les cœurs
    python balayage_grands_R.py --proc 4        # limiter à 4 processus
    python balayage_grands_R.py --sans-super    # ellipses et disque seulement

Le calcul est repris là où il s'est arrêté : chaque paquet de 16 tirages terminé
est ajouté immédiatement à balayage_grands_R.csv, et les paquets déjà présents
sont sautés au relancement (on peut donc l'interrompre par Ctrl+C).
"""
import sys, os, time, math, zlib, argparse
import numpy as np  # pyright: ignore[reportMissingImports]
from multiprocessing import Pool
from balayage_translations import columns, hull_perim, invariants

N_SAMPLES = 128
R_VALUES = [30000, 46000, 70000, 110000, 170000, 260000, 400000, 630000, 1000000]
SHAPES = [("ellipse", 1.0), ("ellipse", 2.0), ("ellipse", 5.0), ("super", 3.0)]
SORTIE = "balayage_grands_R.csv"
ENTETE = "forme,param,R,echantillon,theta,u,v,P_cont,P_hull,deficit_abs,m_sommets\n"


PAQUET = 16  # tirages par tâche, pour répartir les gros R sur tous les cœurs


def tirages(kind, prm, R):
    """Les N_SAMPLES tirages (θ, u, v) d'un point, toujours les mêmes."""
    rng = np.random.default_rng(zlib.crc32(f"grandR-{kind}-{prm}-{int(R)}".encode()))
    return [(rng.uniform(0, math.pi / 2), *rng.uniform(0, 1, 2)) for _ in range(N_SAMPLES)]


def job(args):
    kind, prm, R, P1, i0 = args
    tir = tirages(kind, prm, R)
    out = []
    for i in range(i0, min(i0 + PAQUET, N_SAMPLES)):
        th, u, v = tir[i]
        Ph, m = hull_perim(columns(kind, prm, float(R), th, u, v))
        Pc = P1 * R
        assert Ph <= Pc + 1e-6, (kind, prm, R, i)
        out.append((kind, prm, R, i, th, u, v, Pc, Ph, Pc - Ph, m))
    return out


def deja_faits():
    """Ensemble des (forme, param, R, échantillon) déjà écrits."""
    if not os.path.exists(SORTIE):
        return set()
    faits = set()
    with open(SORTIE) as f:
        next(f)
        for ligne in f:
            c = ligne.split(",")
            if len(c) == 11:
                faits.add((c[0], float(c[1]), int(float(c[2])), int(c[3])))
    return faits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proc", type=int, default=None)
    ap.add_argument("--sans-super", action="store_true")
    a = ap.parse_args()
    shapes = [s for s in SHAPES if not (a.sans_super and s[0] == "super")]
    inv = {sh: invariants(*sh) for sh in shapes}
    faits = deja_faits()
    tasks = [(k, p, R, inv[(k, p)][0], i0) for (k, p) in shapes for R in R_VALUES
             for i0 in range(0, N_SAMPLES, PAQUET)
             if not all((k, p, R, i) in faits
                        for i in range(i0, min(i0 + PAQUET, N_SAMPLES)))]
    # les plus coûteux d'abord (superellipse, grands R)
    tasks.sort(key=lambda t: -(t[2] * (10 if t[0] == "super" else 1)))
    print(f"{len(tasks)} paquets de {PAQUET} tirages à calculer"
          f" ({len(faits)} tirages déjà faits).", flush=True)
    if not os.path.exists(SORTIE):
        with open(SORTIE, "w") as f:
            f.write(ENTETE)
    t0 = time.time()
    with Pool(a.proc) as pool:
        for n, res in enumerate(pool.imap_unordered(job, tasks), 1):
            with open(SORTIE, "a") as f:
                for r in res:
                    f.write(",".join(str(x) for x in r) + "\n")
            k, p, R, i = res[0][0], res[0][1], res[0][2], res[0][3]
            print(f"  {n}/{len(tasks)}  {k} {p:g}  R = {R:>7d}  tirages {i}-{i+len(res)-1}"
                  f"  ({time.time()-t0:.0f} s)", flush=True)
    print(f"Terminé en {time.time()-t0:.0f} s. Lancer ensuite : python analyse_grands_R.py")


if __name__ == "__main__":
    main()
