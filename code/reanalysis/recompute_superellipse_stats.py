import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_EXTERNAL = Path(
    r"C:\Users\dynagenese\Desktop\Extension à d'autres constantes transcendantes irrationnelles"
) / "scripts de développements"
DEFAULT_CANDIDATE = SCRIPT_DIR.parents[1] / "data" / "raw" if len(SCRIPT_DIR.parents) >= 2 else DEFAULT_EXTERNAL
DATA_DIR = DEFAULT_CANDIDATE if DEFAULT_CANDIDATE.exists() else DEFAULT_EXTERNAL
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(DATA_DIR))

import enveloppe_commune as EC  # noqa: E402


def gamma_pred(p: float) -> float:
    return 2 * (4 * p - 5) / (3 * (p + 1))


def recompute(p: float, data_dir: Path = DATA_DIR) -> None:
    source = data_dir / f"superellipse_p{p}_singuliere.json"
    raw = json.loads(source.read_text(encoding="utf-8"))
    perimeter, i43, i32 = EC.invariants_superellipse(p)
    ratio = i32 / i43
    rows = []
    for radius in sorted(map(int, raw)):
        values = np.array(list(raw[str(radius)].values()), dtype=float)
        factor = radius ** (1 / 3) / i43
        c_eff = (perimeter * radius - values.mean()) * factor
        stderr = values.std(ddof=1) / math.sqrt(len(values)) * factor
        c_theory = EC.C_theo(radius, ratio)
        rows.append((radius, len(values), c_eff, stderr, c_theory))

    radii = np.array([r[0] for r in rows], dtype=float)
    delta = np.array([r[2] - r[4] for r in rows])
    sigma = np.array([r[3] for r in rows])
    gamma0 = gamma_pred(p)
    weights = 1 / sigma**2
    x0 = radii ** (-gamma0)
    amplitude0 = float(np.sum(weights * delta * x0) / np.sum(weights * x0**2))
    amp0_error = float(1 / math.sqrt(np.sum(weights * x0**2)))
    chi_null = float(np.sum((delta / sigma) ** 2))
    chi_fixed = float(np.sum(((delta - amplitude0 * x0) / sigma) ** 2))
    (amplitude1, gamma1), covariance = curve_fit(
        lambda x, amplitude, gamma: amplitude * x ** (-gamma),
        radii,
        delta,
        p0=[amplitude0, gamma0],
        sigma=sigma,
        absolute_sigma=True,
        maxfev=20000,
    )
    amp1_error, gamma1_error = np.sqrt(np.diag(covariance))
    chi_free = float(
        np.sum(((delta - amplitude1 * radii ** (-gamma1)) / sigma) ** 2)
    )

    print(f"source={source}")
    print(f"p={p}; P={perimeter:.10f}; I43={i43:.10f}; I32={i32:.10f}")
    print("R,n,Ceff,SE,Ctheory,residual_sigma")
    for radius, n, c_eff, stderr, c_theory in rows:
        print(
            f"{radius},{n},{c_eff:.8f},{stderr:.8f},{c_theory:.8f},"
            f"{(c_eff-c_theory)/stderr:+.4f}"
        )
    print(f"chi_null={chi_null:.6f}; n={len(rows)}")
    print(
        f"fixed_gamma={gamma0:.8f}; A={amplitude0:.8f}+/-{amp0_error:.8f}; "
        f"chi={chi_fixed:.6f}; dof={len(rows)-1}"
    )
    print(
        f"free_gamma={gamma1:.8f}+/-{gamma1_error:.8f}; "
        f"A={amplitude1:.8f}+/-{amp1_error:.8f}; "
        f"chi={chi_free:.6f}; dof={len(rows)-2}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=float, nargs="?", default=1.75)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    options = parser.parse_args()
    recompute(options.p, options.data_dir)
