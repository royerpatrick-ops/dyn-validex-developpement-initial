"""Reanalyse locale des sorties existantes pour la candidate Zenodo.

Ce script ne genere aucune simulation. Il lit uniquement les CSV historiques,
regroupe les repetitions par (forme, parametre, R), recalcule les erreurs-types
et teste la loi a deux termes sur les 92 points annonces (R >= 1230).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad


TARGETS = {("ellipse", 1.0), ("ellipse", 2.0), ("ellipse", 5.0), ("super", 3.0)}
C_INF = 0.7192
B_CONST = 0.2894


def i_super(p: float, exponent: float) -> float:
    def f(x: float) -> float:
        if x <= 0.0:
            return 0.0 if exponent > 0 else 1.0
        w = 1.0 - x**p
        yp = -(x ** (p - 1)) * (w ** (1 / p - 1))
        ypp = -(p - 1) * (x ** (p - 2)) * (w ** (1 / p - 2))
        g = 1.0 + yp * yp
        return (abs(ypp) / g**1.5) ** exponent * math.sqrt(g)

    xm = 2.0 ** (-1.0 / p)
    value, _ = quad(
        lambda s: f(s**6) * 6 * s**5,
        0.0,
        xm ** (1 / 6),
        limit=2000,
        epsabs=1e-13,
        epsrel=1e-12,
    )
    return 8.0 * value


def i_ellipse(k: float, exponent: float) -> float:
    major, minor = 1.0, 1.0 / k

    def speed(t: float) -> float:
        return math.hypot(major * math.sin(t), minor * math.cos(t))

    value, _ = quad(
        lambda t: (major * minor / speed(t) ** 3) ** exponent * speed(t),
        0.0,
        math.pi / 2,
        limit=1000,
        epsabs=1e-13,
        epsrel=1e-12,
    )
    return 4.0 * value


def invariants(kind: str, param: float) -> tuple[float, float]:
    fn = i_super if kind == "super" else i_ellipse
    return fn(param, 4 / 3), fn(param, 3 / 2)


def weighted_fit(y: np.ndarray, se: np.ndarray, design: np.ndarray) -> dict:
    weighted_design = design / se[:, None]
    weighted_y = y / se
    beta, *_ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
    covariance = np.linalg.inv(weighted_design.T @ weighted_design)
    residual = (y - design @ beta) / se
    chi2 = float(np.sum(residual**2))
    return {
        "beta": beta.tolist(),
        "standard_errors_nominal": np.sqrt(np.diag(covariance)).tolist(),
        "covariance_nominal": covariance.tolist(),
        "chi2": chi2,
        "n": int(len(y)),
        "dof": int(len(y) - design.shape[1]),
        "chi2_per_dof": chi2 / (len(y) - design.shape[1]),
    }


def summarize_subset(frame: pd.DataFrame) -> dict:
    """Return descriptive fixed-law and diagonal-WLS statistics for a subset."""
    y = frame.C_eff.to_numpy()
    se = frame.C_eff_se.to_numpy()
    x = (frame.ratio_I32_I43 * frame.R ** (-1 / 6)).to_numpy()
    fixed_chi2 = float(np.sum(frame.z_fixed**2))
    free = weighted_fit(y, se, np.column_stack([np.ones(len(y)), -x]))
    free["parameter_order"] = ["C_inf", "B"]
    return {
        "n": int(len(frame)),
        "fixed_chi2": fixed_chi2,
        "fixed_chi2_per_n": fixed_chi2 / len(frame),
        "fixed_max_abs_z": float(np.max(np.abs(frame.z_fixed))),
        "two_parameter_weighted_fit_diagonal": free,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (root / "data" / "raw" / "balayage_translations.csv").exists():
        data_dir = root / "data" / "raw"
    else:
        data_dir = root / "programmes de developpements"
        if not data_dir.exists():
            data_dir = root / "programmes de développements"

    d1 = pd.read_csv(data_dir / "balayage_translations.csv")
    d2 = pd.read_csv(data_dir / "balayage_grands_R.csv")
    disk = (d1["forme"] == "super") & (d1["param"] == 2.0)
    d1.loc[disk, ["forme", "param"]] = ["ellipse", 1.0]
    raw = pd.concat([d1, d2], ignore_index=True)
    target_mask = [
        (str(kind), float(param)) in TARGETS
        for kind, param in zip(raw["forme"], raw["param"])
    ]
    raw = raw.loc[target_mask].copy()
    raw = raw[raw["R"] >= 1230].copy()

    inv = {key: invariants(*key) for key in TARGETS}
    raw["I43"] = [inv[(str(f), float(p))][0] for f, p in zip(raw.forme, raw.param)]
    raw["I32"] = [inv[(str(f), float(p))][1] for f, p in zip(raw.forme, raw.param)]
    raw["C_eff"] = raw.deficit_abs * raw.R ** (1 / 3) / raw.I43

    grouped = (
        raw.groupby(["forme", "param", "R"], as_index=False)
        .agg(
            C_eff=("C_eff", "mean"),
            C_eff_sd=("C_eff", "std"),
            n=("C_eff", "size"),
            I43=("I43", "first"),
            I32=("I32", "first"),
        )
        .sort_values(["forme", "param", "R"])
    )
    grouped["C_eff_se"] = grouped.C_eff_sd / np.sqrt(grouped.n)
    grouped["ratio_I32_I43"] = grouped.I32 / grouped.I43
    grouped["C_pred_fixed"] = C_INF - B_CONST * grouped.ratio_I32_I43 * grouped.R ** (-1 / 6)
    grouped["z_fixed"] = (grouped.C_eff - grouped.C_pred_fixed) / grouped.C_eff_se

    fixed_chi2 = float(np.sum(grouped.z_fixed**2))
    fixed = {
        "C_inf": C_INF,
        "B": B_CONST,
        "chi2": fixed_chi2,
        "n": int(len(grouped)),
        "chi2_per_n": fixed_chi2 / len(grouped),
        "dof_if_external_constants": int(len(grouped)),
        "max_abs_z": float(np.max(np.abs(grouped.z_fixed))),
        "covariance_assumption": "diagonal; cross-radius covariance not estimable from supplied outputs",
    }

    y = grouped.C_eff.to_numpy()
    se = grouped.C_eff_se.to_numpy()
    x = (grouped.ratio_I32_I43 * grouped.R ** (-1 / 6)).to_numpy()
    free = weighted_fit(y, se, np.column_stack([np.ones(len(y)), -x]))
    free["parameter_order"] = ["C_inf", "B"]
    fixed_c = weighted_fit(y - C_INF, se, (-x)[:, None])
    fixed_c["parameter_order"] = ["B_with_C_inf_fixed"]

    exploratory_mask = (grouped.forme == "super") & (grouped.param == 3.0)
    regular_core = summarize_subset(grouped.loc[~exploratory_mask])
    exploratory_p3 = summarize_subset(grouped.loc[exploratory_mask])
    by_shape = {}
    for (kind, param), frame in grouped.groupby(["forme", "param"], sort=True):
        by_shape[f"{kind}:{param:g}"] = summarize_subset(frame)

    counts = (
        grouped.groupby(["forme", "param"], as_index=False)
        .agg(points=("R", "size"), R_min=("R", "min"), R_max=("R", "max"), n_min=("n", "min"), n_max=("n", "max"))
        .to_dict(orient="records")
    )
    result = {
        "sources": [
            str(data_dir / "balayage_translations.csv"),
            str(data_dir / "balayage_grands_R.csv"),
        ],
        "selection": "four shapes present in grand-R campaign; grouped points with R >= 1230; duplicate R=30000 pooled",
        "uncertainty": "plus/minus is one standard error of the sample mean: sample SD with ddof=1 divided by sqrt(n)",
        "counts_by_shape": counts,
        "fixed_external_constants": fixed,
        "two_parameter_weighted_fit_diagonal": free,
        "one_parameter_weighted_fit_diagonal": fixed_c,
        "class_stratification": {
            "regular_positive_curvature_core": {
                "definition": "disk and ellipses 2:1 and 5:1; strictly positive curvature",
                **regular_core,
            },
            "exploratory_superellipse_p3": {
                "definition": "superellipse p=3; curvature vanishes at the four axial points",
                **exploratory_p3,
            },
            "by_shape": by_shape,
        },
        "notes": [
            "No large generator was rerun.",
            "Nominal covariance matrices are inverse weighted normal matrices and do not include model discrepancy.",
            "The 92 grouped means are not 92 raw draws; n is 64, 128, or 192 at pooled R=30000.",
            "Model and cutoff were selected in the same research workflow; chi-square is descriptive, not a preregistered confirmatory p-value.",
            "The superellipse p=3 is retained in the 92-point aggregate only as an exploratory extension outside the positive-curvature core.",
        ],
    }

    grouped.to_csv(output / "reanalysis_92_points.csv", index=False)
    (output / "reanalysis_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        "# Reanalyse des 92 points",
        "",
        "Aucune campagne numerique n'a ete relancee. Les deux CSV historiques ont ete groupes par forme, parametre et rayon.",
        "",
        f"- Loi fixee C_inf={C_INF}, B={B_CONST}: chi2={fixed_chi2:.6f}; n={len(grouped)}; chi2/n={fixed_chi2/len(grouped):.6f}; max |z|={fixed['max_abs_z']:.4f}.",
        f"- Ajustement libre diagonal: C_inf={free['beta'][0]:.8f} +/- {free['standard_errors_nominal'][0]:.8f}; B={free['beta'][1]:.8f} +/- {free['standard_errors_nominal'][1]:.8f}; chi2/dof={free['chi2_per_dof']:.6f}.",
        f"- C_inf fixe: B={fixed_c['beta'][0]:.8f} +/- {fixed_c['standard_errors_nominal'][0]:.8f}; chi2/dof={fixed_c['chi2_per_dof']:.6f}.",
        f"- Noyau a courbure strictement positive (disque et ellipses 2:1, 5:1): n={regular_core['n']}; chi2 fixe={regular_core['fixed_chi2']:.6f}; chi2/n={regular_core['fixed_chi2_per_n']:.6f}.",
        f"- Extension exploratoire superellipse p=3 (courbure nulle aux axes): n={exploratory_p3['n']}; chi2 fixe={exploratory_p3['fixed_chi2']:.6f}; chi2/n={exploratory_p3['fixed_chi2_per_n']:.6f}.",
        "- Ici +/- des points = erreur-type de moyenne, SD(ddof=1)/sqrt(n).",
        "- La covariance inter-rayons n'est pas identifiable dans les sorties fournies; le calcul utilise donc une covariance diagonale.",
        "- La selection du modele et du seuil R>=1230 n'etant pas preregistree, ces chi2 sont descriptifs.",
        "",
    ]
    (output / "REANALYSE_92_POINTS.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
