"""
Candidate v0.3 : déficit de périmètre de l'enveloppe entière avec moyenne
sur rotations ET translations aléatoires du centre.

Correction par rapport à balayage_corrige.py : pour le disque (et les formes
presque circulaires), une rotation autour d'un centre entier ne change pas
l'ensemble des points du réseau ; la moyenne sur les rotations seules ne
moyenne donc rien. On tire ici, pour chaque échantillon, une rotation
θ ∈ [0, π/2) ET un décalage (u, v) ∈ [0, 1)² du réseau.

Sorties : balayage_translations.csv
  forme, param, R, echantillon, theta, u, v, P_cont, P_hull, deficit_abs, m_sommets
"""
import sys, time, math, zlib
import numpy as np  # pyright: ignore[reportMissingImports]
from multiprocessing import Pool

N_SAMPLES = 64
R_VALUES = np.unique(np.round(np.geomspace(100, 30000, 26)).astype(int)).astype(float)
SHAPES = ([("super", p) for p in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]]
          + [("ellipse", k) for k in [1.5, 2.0, 3.0, 5.0]])


# ------------------------------------------------------------------
# Invariants numeriques (domaine de finitude controle) I_a = ∫ κ^a ds  (superellipse : paramétrage par x,
# changement x = s^6 qui supprime les singularités de courbure)
# ------------------------------------------------------------------
from invariant_checks import I_super, I_ellipse, invariants


# ------------------------------------------------------------------
# Points extrêmes de chaque colonne du réseau Z² - (u, v) dans R·K
# ------------------------------------------------------------------
from geometry_columns import columns, column_indices, integer_hull_perimeter


def hull_perim(pts):
    """v0.3: integer indices required; translated float coordinates are display-only."""
    pts = np.asarray(pts)
    if not np.issubdtype(pts.dtype, np.integer):
        raise ValueError("Use column_indices(), not translated float coordinates")
    return integer_hull_perimeter(pts)


# ------------------------------------------------------------------
def job(args):
    kind, prm, R, P1 = args
    rng = np.random.default_rng(zlib.crc32(f"{kind}-{prm}-{int(R)}".encode()))
    out = []
    for i in range(N_SAMPLES):
        th = rng.uniform(0, math.pi / 2); u, v = rng.uniform(0, 1, 2)
        Ph, m = integer_hull_perimeter(column_indices(kind, prm, R, th, u, v))
        Pc = P1 * R
        assert Ph <= Pc + 1e-6
        out.append((kind, prm, R, i, th, u, v, Pc, Ph, Pc - Ph, m))
    return out


def main():
    inv = {sh: invariants(*sh) for sh in SHAPES}
    tasks = [(k, p, R, inv[(k, p)][0]) for (k, p) in SHAPES for R in R_VALUES]
    tasks.sort(key=lambda t: -t[2])  # gros R d'abord (équilibrage)
    t0 = time.time(); rows = []
    with Pool() as pool:
        for n, res in enumerate(pool.imap_unordered(job, tasks), 1):
            rows.extend(res)
            if n % 20 == 0:
                print(f"{n}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)
    with open("balayage_translations.csv", "w") as f:
        f.write("forme,param,R,echantillon,theta,u,v,P_cont,P_hull,deficit_abs,m_sommets\n")
        for r in sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3])):
            f.write(",".join(str(x) for x in r) + "\n")
    with open("invariants.csv", "w") as f:
        f.write("forme,param,P,Omega,I43\n")
        for (k, p), (P, Om, J) in inv.items():
            f.write(f"{k},{p},{P:.12f},{Om:.12f},{J:.12f}\n")
    print(f"Terminé en {time.time()-t0:.0f}s")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Campagne originale, invariants corriges")
    parser.add_argument("--campaign", action="store_true", help="lance explicitement la campagne complete")
    parser.add_argument("--out", type=str, help="nouveau dossier pour la campagne")
    args = parser.parse_args()
    if args.campaign:
        parser.error("Campagnes suspendues dans v0.3: audit geometrique requis")
        if not args.out:
            parser.error("--out requis avec --campaign")
        import os
        from pathlib import Path
        output_dir = Path(args.out).resolve()
        output_dir.mkdir(parents=True, exist_ok=False)
        os.chdir(output_dir)
        main()
    else:
        parser.print_help()
