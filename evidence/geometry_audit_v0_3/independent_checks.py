from __future__ import annotations

import ast
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
import time
import traceback
import zipfile
from fractions import Fraction as F
from pathlib import Path, PurePosixPath

import mpmath
import numpy as np
import scipy
from mpmath.ctx_mp import MPContext


ROOT = Path(__file__).resolve().parents[1]
AUDIT = Path(__file__).resolve().parent
CANDIDATE = ROOT / "candidate"
PYTHON = Path(sys.executable).resolve()
sys.path.insert(0, str(CANDIDATE))

from geometry_columns import (  # noqa: E402
    column_indices,
    disk_column_indices,
    integer_hull_perimeter,
    integer_hull_vertices,
)
from geometry_predicates import (  # noqa: E402
    BodyPredicate,
    GeometryIndeterminate,
    endpoints,
    rational,
)


OUTPUT = AUDIT / "independent_results.json"
if OUTPUT.exists():
    raise SystemExit(f"refus d'ecraser {OUTPUT}")

started = time.monotonic()
checks = []
observations = []


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def floor_q(value: F) -> int:
    return value.numerator // value.denominator


def ceil_q(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def run_check(identifier: str, title: str, function) -> None:
    begin = time.monotonic()
    try:
        details = function()
        checks.append(
            {
                "id": identifier,
                "title": title,
                "status": "PASS",
                "seconds": time.monotonic() - begin,
                "details": details,
            }
        )
    except Exception as exc:
        checks.append(
            {
                "id": identifier,
                "title": title,
                "status": "FAIL",
                "seconds": time.monotonic() - begin,
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
        )


def extrema(points) -> set[tuple[int, int]]:
    by_x = {}
    for x, y in points:
        x, y = int(x), int(y)
        if x not in by_x:
            by_x[x] = [y, y]
        else:
            by_x[x][0] = min(by_x[x][0], y)
            by_x[x][1] = max(by_x[x][1], y)
    return {(x, y) for x, pair in by_x.items() for y in set(pair)}


def jarvis(points) -> set[tuple[int, int]]:
    unique = sorted(set((int(x), int(y)) for x, y in points))
    if len(unique) < 3:
        return set(unique)
    start = unique[0]
    current = start
    result = []
    for _ in range(len(unique) + 1):
        result.append(current)
        candidate = next(point for point in unique if point != current)
        for point in unique:
            if point == current:
                continue
            cross = (
                (candidate[0] - current[0]) * (point[1] - current[1])
                - (candidate[1] - current[1]) * (point[0] - current[0])
            )
            candidate_distance = (candidate[0] - current[0]) ** 2 + (
                candidate[1] - current[1]
            ) ** 2
            point_distance = (point[0] - current[0]) ** 2 + (
                point[1] - current[1]
            ) ** 2
            if cross < 0 or (cross == 0 and point_distance > candidate_distance):
                candidate = point
        current = candidate
        if current == start:
            return set(result)
    raise AssertionError("Jarvis n'a pas referme l'enveloppe")


def brute_disk(radius, u=0, v=0):
    radius, u, v = F(radius), F(u), F(v)
    points = []
    for x in range(ceil_q(u - radius), floor_q(u + radius) + 1):
        for y in range(ceil_q(v - radius), floor_q(v + radius) + 1):
            if (F(x) - u) ** 2 + (F(y) - v) ** 2 <= radius**2:
                points.append((x, y))
    return points


def mpq(context: MPContext, value) -> object:
    value = F(value)
    return context.mpf(value.numerator) / value.denominator


def exact_axis_value(kind, parameter: F, radius: F, x: F, y: F):
    X, Y = abs(x / radius), abs(y / radius)
    if kind == "ellipse":
        return X * X + (parameter * Y) ** 2
    if parameter.denominator == 1:
        return X ** int(parameter) + Y ** int(parameter)
    return None


def high_value(context, kind, parameter, radius, theta, u, v, x, y):
    parameter, radius, theta, u, v = map(F, (parameter, radius, theta, u, v))
    x, y = F(x) - u, F(y) - v
    if (kind == "ellipse" and parameter == 1) or (
        kind == "super" and parameter == 2
    ):
        return mpq(context, (x * x + y * y) / (radius * radius))
    c, s = context.cos(mpq(context, theta)), context.sin(mpq(context, theta))
    X = abs((c * mpq(context, x) + s * mpq(context, y)) / mpq(context, radius))
    Y = abs((-s * mpq(context, x) + c * mpq(context, y)) / mpq(context, radius))
    if kind == "ellipse":
        return X**2 + (mpq(context, parameter) * Y) ** 2
    power = mpq(context, parameter)
    return X**power + Y**power


def brute_shape(kind, parameter, radius, theta, u, v):
    kind = str(kind)
    parameter, radius, theta, u, v = map(F, (parameter, radius, theta, u, v))
    bound = (
        max(radius, radius / parameter) if kind == "ellipse" else 2 * radius
    )
    context = MPContext()
    context.dps = 220
    threshold = context.mpf("1e-180")
    points = []
    ambiguous = []
    for x in range(ceil_q(u - bound), floor_q(u + bound) + 1):
        for y in range(ceil_q(v - bound), floor_q(v + bound) + 1):
            relative_x, relative_y = F(x) - u, F(y) - v
            exact = None
            if (kind == "ellipse" and parameter == 1) or (
                kind == "super" and parameter == 2
            ):
                exact = (relative_x**2 + relative_y**2) / radius**2
            elif theta == 0:
                exact = exact_axis_value(
                    kind, parameter, radius, relative_x, relative_y
                )
            if exact is not None:
                inside = exact <= 1
            else:
                value = high_value(
                    context, kind, parameter, radius, theta, u, v, x, y
                )
                if abs(value - 1) < threshold:
                    ambiguous.append((x, y, context.nstr(value - 1, 12)))
                    continue
                inside = value < 1
            if inside:
                points.append((x, y))
    require(not ambiguous, f"reference haute precision ambigue: {ambiguous[:3]}")
    return points


def check_integrity():
    package_manifest = json.loads((ROOT / "package_manifest.json").read_text())
    mismatches = []
    for relative, expected in package_manifest.items():
        path = ROOT / relative
        if not path.is_file():
            mismatches.append({"path": relative, "error": "missing"})
        elif sha256(path) != expected:
            mismatches.append(
                {
                    "path": relative,
                    "expected": expected,
                    "actual": sha256(path),
                }
            )
    require(not mismatches, f"manifest package: {mismatches[:3]}")
    originals_manifest = json.loads((ROOT / "originals_sha256.json").read_text())
    originals = {
        str(path.relative_to(ROOT / "originals")).replace("\\", "/")
        for path in (ROOT / "originals").rglob("*")
        if path.is_file()
    }
    require(len(originals_manifest) == 89, "nombre de sources originales inattendu")
    require(originals == set(originals_manifest), "inventaire originals divergent")
    for relative, expected in originals_manifest.items():
        require(
            sha256(ROOT / "originals" / relative) == expected,
            f"original modifie: {relative}",
        )
    expected_candidate = {
        Path(relative).name: digest
        for relative, digest in package_manifest.items()
        if relative.startswith("candidate/")
    }
    actual_candidate = {
        path.name: sha256(path) for path in sorted(CANDIDATE.glob("*.py"))
    }
    require(actual_candidate == expected_candidate, "candidate divergent du manifeste")
    return {
        "package_manifest_entries": len(package_manifest),
        "original_files": len(originals_manifest),
        "candidate_files": len(actual_candidate),
        "candidate_sha256": actual_candidate,
    }


def compare_archive(archive: Path, disk_root: Path):
    mismatches = []
    with zipfile.ZipFile(archive) as handle:
        names = [name for name in handle.namelist() if not name.endswith("/")]
        for name in names:
            parts = PurePosixPath(name).parts
            relative = Path(*parts[1:])
            local = disk_root / relative
            if not local.is_file():
                mismatches.append({"path": str(relative), "error": "missing"})
            else:
                archived_hash = hashlib.sha256(handle.read(name)).hexdigest()
                if sha256(local) != archived_hash:
                    mismatches.append(
                        {
                            "path": str(relative),
                            "archive": archived_hash,
                            "disk": sha256(local),
                        }
                    )
    require(not mismatches, f"archive historique divergente: {mismatches[:3]}")
    return len(names)


def check_history():
    parent = ROOT.parents[1]
    v02 = parent / "DYN-CORRECTIVE_INTEGER-HULL_v0_2" / "DYN-CORRECTIVE_INTEGER-HULL_v0_2"
    review = (
        parent
        / "DYN-REVUE_AUDIT-ET-PATCH_INTEGER-HULL_v0_1"
        / "DYN-REVUE_AUDIT-ET-PATCH_INTEGER-HULL_v0_1"
    )
    require(v02.is_dir(), f"paquet v0.2 absent: {v02}")
    require(review.is_dir(), f"revue v0.1 absente: {review}")
    return {
        "v0_2_files_compared": compare_archive(
            ROOT / "history" / "DYN-CORRECTIVE_INTEGER-HULL_v0_2.zip", v02
        ),
        "review_v0_1_files_compared": compare_archive(
            ROOT / "history" / "DYN-REVUE_AUDIT-ET-PATCH_INTEGER-HULL_v0_1.zip",
            review,
        ),
    }


def check_ast():
    parsed = []
    for path in sorted(CANDIDATE.glob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        parsed.append(path.name)
    return {"parsed": parsed, "count": len(parsed)}


def check_inputs():
    require(rational(0.1) == F(0.1), "float non converti en rationnel binaire exact")
    require(rational(np.float64(0.1)) == F(float(np.float64(0.1))), "numpy float")
    require(rational(np.int64(7)) == 7, "numpy int")
    invalid = [True, False, float("nan"), float("inf"), "1", None]
    rejected = 0
    for value in invalid:
        try:
            rational(value)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError(f"entree invalide acceptee: {value!r}")
    for arguments in [
        ("wrong", 1, 2, 0, 0, 0),
        ("ellipse", 0, 2, 0, 0, 0),
        ("super", 1, 2, 0, 0, 0),
        ("ellipse", 1, 2, 0, 1, 0),
    ]:
        try:
            BodyPredicate(*arguments)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError(f"parametres invalides acceptes: {arguments}")
    return {"invalid_inputs_rejected": rejected, "float_semantics": "exact binary rational"}


def check_exact_disks():
    radii = [F(1, 7), F(2, 3), F(7, 6), F(5, 2), F(17, 3), F(31, 5)]
    shifts = [
        (F(0), F(0)),
        (F(1, 3), F(2, 5)),
        (F(2, 7), F(4, 9)),
        (F(1, 97), F(96, 101)),
        (F(6, 7), F(1, 11)),
    ]
    cases = 0
    lattice_points = 0
    for radius in radii:
        for u, v in shifts:
            reference = brute_disk(radius, u, v)
            actual = disk_column_indices(radius, u, v)
            require(
                set(map(tuple, actual)) == extrema(reference),
                f"extremes disque R={radius},u={u},v={v}",
            )
            actual_hull = set(map(tuple, integer_hull_vertices(actual)))
            require(
                actual_hull == jarvis(reference),
                f"enveloppe disque R={radius},u={u},v={v}",
            )
            cases += 1
            lattice_points += len(reference)
    return {
        "cases": cases,
        "full_grid_points_classified": lattice_points,
        "denominators": [3, 5, 7, 9, 11, 97, 101],
    }


def check_float_disks():
    smallest = float(np.nextafter(0.0, 1.0))
    radii = [
        math.nextafter(3.0, 0.0),
        3.0,
        math.nextafter(3.0, math.inf),
        math.nextafter(5.5, 0.0),
        5.5,
        math.nextafter(5.5, math.inf),
    ]
    shifts = [
        (smallest, 0.0),
        (float(1 / 3), float(2 / 3)),
        (math.nextafter(0.5, 0.0), math.nextafter(0.5, math.inf)),
    ]
    cases = 0
    for radius in radii:
        for u, v in shifts:
            reference = brute_disk(F(radius), F(u), F(v))
            actual = column_indices("ellipse", 1, radius, 0.314159, u, v)
            require(
                set(map(tuple, actual)) == extrema(reference),
                f"disque float R={radius!r},u={u!r},v={v!r}",
            )
            cases += 1
    for radius in [F(31), F(31, 3), F(17, 5)]:
        body = BodyPredicate("ellipse", 1, radius, F(314159, 1000000))
        require(body.contains_point(radius, 0), f"frontiere exclue R={radius}")
        require(
            not body.contains_point(radius + F(1, 10**120), 0),
            f"exterieur immediat inclus R={radius}",
        )
    return {"float_cases": cases, "exact_boundary_pairs": 3, "subnormal_shift": smallest}


def check_general_shapes():
    families = [
        ("super", F(3, 2)),
        ("super", F(5, 2)),
        ("super", F(3)),
        ("super", F(4)),
        ("super", F(5)),
        ("super", F(10)),
        ("ellipse", F(1)),
        ("ellipse", F(3, 2)),
        ("ellipse", F(2)),
        ("ellipse", F(3)),
        ("ellipse", F(5)),
    ]
    configurations = [
        (F(0), F(1, 3), F(2, 5)),
        (F(314159, 1000000), F(2, 7), F(4, 9)),
        (F(-2718, 10000), F(1, 97), F(5, 11)),
    ]
    cases = 0
    classified = 0
    for kind, parameter in families:
        for theta, u, v in configurations:
            radius = F(13, 3)
            reference = brute_shape(kind, parameter, radius, theta, u, v)
            actual = column_indices(kind, parameter, radius, theta, u, v)
            require(
                set(map(tuple, actual)) == extrema(reference),
                f"colonnes {kind}/{parameter},theta={theta},u={u},v={v}",
            )
            require(
                set(map(tuple, integer_hull_vertices(actual))) == jarvis(reference),
                f"enveloppe {kind}/{parameter},theta={theta},u={u},v={v}",
            )
            cases += 1
            classified += len(reference)
    return {
        "families": len(families),
        "configurations_per_family": len(configurations),
        "cases": cases,
        "inside_lattice_points": classified,
    }


def check_near_generic_boundaries():
    context = MPContext()
    context.dps = 240
    cases = [
        ("ellipse", F(3, 2), 0.314159, 1 / 3, 2 / 5, 3, 2),
        ("ellipse", F(5), -0.2718, 2 / 7, 4 / 9, -2, 3),
        ("super", F(3, 2), 0.314159, 1 / 97, 5 / 11, 2, -3),
        ("super", F(13, 10), -0.2718, 2 / 5, 1 / 7, -3, -2),
        ("super", F(5, 2), 0.71, 1 / 3, 2 / 7, 4, 1),
        ("super", F(10), -0.19, 3 / 8, 7 / 9, 1, 4),
    ]
    tested = 0
    refinements = 0
    for kind, parameter, theta, u, v, x, y in cases:
        theta_q, u_q, v_q = F(theta), F(u), F(v)
        c, s = context.cos(mpq(context, theta_q)), context.sin(mpq(context, theta_q))
        relative_x, relative_y = mpq(context, F(x) - u_q), mpq(context, F(y) - v_q)
        X = abs(c * relative_x + s * relative_y)
        Y = abs(-s * relative_x + c * relative_y)
        if kind == "ellipse":
            critical = context.sqrt(X**2 + (mpq(context, parameter) * Y) ** 2)
        else:
            power = mpq(context, parameter)
            critical = (X**power + Y**power) ** (1 / power)
        center = float(critical)
        for radius in [
            math.nextafter(center, 0.0),
            center,
            math.nextafter(center, math.inf),
        ]:
            value = high_value(
                context, kind, parameter, F(radius), theta_q, u_q, v_q, x, y
            )
            require(abs(value - 1) > context.mpf("1e-190"), "reference trop proche")
            expected = value < 1
            body = BodyPredicate(kind, parameter, radius, theta_q, u_q, v_q)
            actual = bool(body.contains_indices(np.array([x]), np.array([y]))[0])
            require(
                actual == expected,
                f"frontiere generique {kind}/{parameter},R={radius!r},point={(x,y)}",
            )
            refinements += body.stats["interval_refinements"]
            tested += 1
    return {"classifications": tested, "interval_refinements": refinements}


def check_float64_enclosures():
    generator = np.random.default_rng(20261002)
    context = MPContext()
    context.dps = 260
    shapes = [
        ("ellipse", F(3, 2)),
        ("ellipse", F(5)),
        ("super", F(3, 2)),
        ("super", F(5, 2)),
        ("super", F(3)),
        ("super", F(10)),
    ]
    tested = 0
    finite = 0
    for kind, parameter in shapes:
        for theta in [F(314159, 1000000), F(-2718, 10000), F(71, 100)]:
            body = BodyPredicate(
                kind,
                parameter,
                F(31, 3),
                theta,
                F(2, 7),
                F(4, 9),
            )
            x = generator.integers(-18, 19, 96, dtype=np.int64)
            y = generator.integers(-18, 19, 96, dtype=np.int64)
            bounds = body.bounds(x, y)
            for index in range(len(x)):
                value = high_value(
                    context,
                    kind,
                    parameter,
                    F(31, 3),
                    theta,
                    F(2, 7),
                    F(4, 9),
                    int(x[index]),
                    int(y[index]),
                )
                low, high = float(bounds.lo[index]), float(bounds.hi[index])
                low_mp = context.ninf if math.isinf(low) and low < 0 else mpq(context, F(low))
                high_mp = context.inf if math.isinf(high) and high > 0 else mpq(context, F(high))
                require(
                    low_mp <= value <= high_mp,
                    f"filtre non englobant {kind}/{parameter}: {low}, {value}, {high}",
                )
                finite += int(math.isfinite(low) and math.isfinite(high))
                tested += 1
    return {"values": tested, "finite_intervals": finite, "reference_digits": 260}


def check_mpmath_intervals():
    generator = np.random.default_rng(19071984)
    reference = MPContext()
    reference.dps = 260
    cases = [
        ("ellipse", F(3, 2)),
        ("super", F(13, 10)),
        ("super", F(3, 2)),
        ("super", F(5, 2)),
        ("super", F(10)),
    ]
    tested = 0
    for kind, parameter in cases:
        body = BodyPredicate(
            kind,
            parameter,
            F(31, 3),
            F(314159, 1000000),
            F(1, 3),
            F(2, 5),
        )
        for _ in range(80):
            x = F(int(generator.integers(-18, 19))) - F(1, 3)
            y = F(int(generator.integers(-18, 19))) - F(2, 5)
            low, high = endpoints(body._iv_value(x, y, 100))
            value = high_value(
                reference,
                kind,
                parameter,
                F(31, 3),
                F(314159, 1000000),
                0,
                0,
                x,
                y,
            )
            require(
                mpq(reference, low) <= value <= mpq(reference, high),
                f"mpmath.iv non englobant {kind}/{parameter}",
            )
            tested += 1
    return {"values": tested, "interval_digits": 100, "reference_digits": 260}


def check_indeterminate():
    context = MPContext()
    context.dps = 100
    theta = context.mpf(1) / 10
    x, y = F(str(context.cos(theta))), F(str(context.sin(theta)))
    body = BodyPredicate("ellipse", 2, 1, F(1, 10), precisions=(15,))
    try:
        body.contains_point(x, y)
    except GeometryIndeterminate:
        pass
    else:
        raise AssertionError("cas ambigu accepte ou rejete silencieusement")
    require(body.stats["indeterminate"] == 1, "compteur d'indetermination")
    return {"natural_indeterminate_cases": 1, "precision_digits": 15}


def check_symmetry_resources_hull():
    points = column_indices("super", F(3, 2), 7, 0, F(1, 2), F(1, 2))
    point_set = set(map(tuple, points))
    reflected = {(1 - int(x), 1 - int(y)) for x, y in point_set}
    require(point_set == reflected, "symetrie centrale rompue")
    try:
        column_indices("ellipse", 1, 100, 0, max_columns=2)
    except ValueError:
        pass
    else:
        raise AssertionError("budget de colonnes non applique")
    try:
        column_indices("ellipse", F(1, 2**50), 1, 0)
    except ValueError:
        pass
    else:
        raise AssertionError("borne d'indices non appliquee")
    square = np.array(
        [
            [-(2**40), -(2**40)],
            [2**40, -(2**40)],
            [2**40, 2**40],
            [-(2**40), 2**40],
            [0, 0],
        ],
        dtype=np.int64,
    )
    vertices = set(map(tuple, integer_hull_vertices(square)))
    require(len(vertices) == 4, "orientation entiere incorrecte")
    perimeter, count = integer_hull_perimeter(square)
    require(perimeter == float(8 * 2**40) and count == 4, "perimetre carre")
    try:
        integer_hull_vertices(np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]))
    except ValueError:
        pass
    else:
        raise AssertionError("indices float acceptes")
    try:
        integer_hull_perimeter(np.array([[0, 0], [1, 1]], dtype=np.int64))
    except ValueError:
        pass
    else:
        raise AssertionError("enveloppe degeneree acceptee")
    bool_budget_accepted = True
    try:
        column_indices("ellipse", 1, F(1, 10), 0, max_columns=True)
    except ValueError:
        bool_budget_accepted = False
    empty_float_accepted = True
    try:
        integer_hull_vertices(np.empty((0, 2), dtype=float))
    except ValueError:
        empty_float_accepted = False
    observations.extend(
        [
            {
                "severity": "LOW",
                "subject": "budget boolean",
                "observed": "accepted" if bool_budget_accepted else "rejected",
            },
            {
                "severity": "LOW",
                "subject": "empty float hull input",
                "observed": "accepted" if empty_float_accepted else "rejected",
            },
        ]
    )
    return {
        "central_symmetry_points": len(point_set),
        "column_budget_rejected": True,
        "coordinate_budget_rejected": True,
        "large_integer_hull_vertices": count,
        "nonempty_float_indices_rejected": True,
        "degenerate_perimeter_rejected": True,
    }


def check_official_outputs():
    unit_status = json.loads((AUDIT / "unittest_status.json").read_text(encoding="utf-8-sig"))
    unit_text = (AUDIT / "unittest_stderr.txt").read_text(encoding="utf-8")
    geometry_status = json.loads((AUDIT / "geometry_status.json").read_text(encoding="utf-8-sig"))
    geometry = json.loads((AUDIT / "geometry" / "result.json").read_text())
    require(unit_status["exit_code"] == 0, "unittest non nul")
    require("Ran 33 tests" in unit_text and unit_text.rstrip().endswith("OK"), "33 tests non prouves")
    require(geometry_status["exit_code"] == 0, "geometry check non nul")
    require(geometry["status"] == "PASS_LOCAL_GEOMETRY_REGRESSIONS", "statut geometry")
    require(geometry["all_radii_2_250"] == {"cases": 1494, "failures": [], "passed": 1494}, "G05")
    shifted = geometry["failed_radii_and_1000_with_shifts"]
    require(shifted == {"cases": 780, "passed": 780, "failures": []}, "G06")
    require(len(geometry["old_csv_replays"]) == 25, "replays historiques")
    require(geometry["old_csv_replay_same_vertex_count"], "comptes sommets historiques")
    require(not geometry["C_recalculated"] and not geometry["full_campaign_run"], "calcul interdit")
    return {
        "unit_tests": 33,
        "all_radii": "1494/1494",
        "shifted_cases": "780/780",
        "historical_replays": 25,
        "max_abs_perimeter_difference": geometry[
            "old_csv_replay_max_abs_perimeter_difference"
        ],
    }


def check_suspensions_and_outputs():
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    geometry_result = AUDIT / "geometry" / "result.json"
    before = sha256(geometry_result)
    refusal = subprocess.run(
        [
            str(PYTHON),
            str(CANDIDATE / "run_geometry_checks.py"),
            "--out",
            str(AUDIT / "geometry"),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        env=environment,
        timeout=60,
    )
    require(refusal.returncode != 0, "sortie existante acceptee")
    require(sha256(geometry_result) == before, "sortie existante modifiee")
    cli_results = []
    with tempfile.TemporaryDirectory(dir=AUDIT) as temporary:
        temporary_path = Path(temporary)
        invocations = [
            (
                "balayage_translations.py",
                ["--campaign", "--out", str(temporary_path / "campaign")],
                temporary_path / "campaign",
            ),
            (
                "formule_integrale.py",
                ["--profile", "reference", "--out", str(temporary_path / "C")],
                temporary_path / "C",
            ),
            (
                "test_classes_w0.py",
                ["--out", str(temporary_path / "classes.json")],
                temporary_path / "classes.json",
            ),
        ]
        for script, arguments, target in invocations:
            completed = subprocess.run(
                [str(PYTHON), str(CANDIDATE / script), *arguments],
                cwd=temporary_path,
                capture_output=True,
                text=True,
                env=environment,
                timeout=60,
            )
            require(completed.returncode != 0, f"CLI non suspendue: {script}")
            require("suspendu" in completed.stderr.lower(), f"message absent: {script}")
            require(not target.exists(), f"sortie partielle: {target}")
            cli_results.append({"script": script, "returncode": completed.returncode})
        import_dir = temporary_path / "imports"
        import_dir.mkdir()
        import_environment = environment.copy()
        import_environment["PYTHONPATH"] = str(CANDIDATE)
        completed = subprocess.run(
            [
                str(PYTHON),
                "-c",
                "import balayage_translations, formule_integrale, test_classes_w0, check_datasets",
            ],
            cwd=import_dir,
            capture_output=True,
            text=True,
            env=import_environment,
            timeout=60,
        )
        require(completed.returncode == 0, f"imports: {completed.stderr}")
        require(not any(import_dir.iterdir()), "effet de bord a l'import")
    return {
        "existing_output_refused": True,
        "existing_output_unchanged": True,
        "suspended_clis": cli_results,
        "imports_without_output": True,
    }


def function_dump(path: Path, names: set[str]):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {
        node.name: ast.dump(node, include_attributes=False)
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names
    }


def check_non_geometry_regressions():
    old = (
        ROOT.parents[1]
        / "DYN-CORRECTIVE_INTEGER-HULL_v0_2"
        / "DYN-CORRECTIVE_INTEGER-HULL_v0_2"
        / "candidate"
    )
    unchanged = [
        "check_datasets.py",
        "coefficient_B.py",
        "integral_model.py",
        "invariant_checks.py",
    ]
    for name in unchanged:
        require(sha256(CANDIDATE / name) == sha256(old / name), f"module modifie: {name}")
    numerical_functions = {"log_nodes", "deterministic_tail", "estimate"}
    require(
        function_dump(CANDIDATE / "formule_integrale.py", numerical_functions)
        == function_dump(old / "formule_integrale.py", numerical_functions),
        "noyau formule_integrale modifie",
    )
    from formule_integrale import G0, SIEGEL
    from test_classes_w0 import K

    require(G0 == 4 * math.sqrt(2) / 15, "G0 modifie")
    require(K == 4 * math.sqrt(2) / 15, "K modifie")
    require(SIEGEL == 6 / math.pi**2, "facteur de Siegel modifie")
    return {
        "byte_identical_modules": unchanged,
        "unchanged_numerical_functions": sorted(numerical_functions),
        "C_calculated": False,
        "campaign_run": False,
    }


run_check("G01a", "Integrite du paquet, des originaux et de la candidate", check_integrity)
run_check("G01b", "Identite des deux archives historiques", check_history)
run_check("G01c", "Analyse syntaxique des douze modules candidats", check_ast)
run_check("G02", "Conversions rationnelles et entrees invalides", check_inputs)
run_check("G03-G04-G08", "Disques rationnels exhaustifs et frontiere", check_exact_disks)
run_check("G03-G04", "Flottants binaires, ULP et subnormal", check_float_disks)
run_check("G09-G10", "Onze familles, toutes les extremites de colonnes", check_general_shapes)
run_check("G04-G09-G12", "Flottants generiques proches de la frontiere", check_near_generic_boundaries)
run_check("G11a", "Encadrements float64 contre 260 chiffres", check_float64_enclosures)
run_check("G11b", "Intervalles mpmath contre 260 chiffres", check_mpmath_intervals)
run_check("G12", "Arret naturel sur cas indetermine", check_indeterminate)
run_check("G13-G14", "Symetries, ressources et enveloppe entiere", check_symmetry_resources_hull)
run_check("G05-G06-G07-G15", "Resultats prescrits et replays historiques", check_official_outputs)
run_check("G16", "Sorties neuves, imports et CLI suspendues", check_suspensions_and_outputs)
run_check("G15", "Noyaux non geometriques preserves", check_non_geometry_regressions)

result = {
    "status": "PASS" if all(check["status"] == "PASS" for check in checks) else "FAIL",
    "scope": "Audit geometrique et regressif fini; aucun calcul de C ni campagne Monte Carlo.",
    "environment": {
        "executable": str(PYTHON),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
    },
    "checks": checks,
    "observations": observations,
    "seconds_total": time.monotonic() - started,
    "C_calculated": False,
    "full_campaign_run": False,
}
with OUTPUT.open("x", encoding="utf-8") as handle:
    json.dump(result, handle, indent=2, ensure_ascii=False)
    handle.write("\n")
print(json.dumps(result, indent=2, ensure_ascii=False))
raise SystemExit(0 if result["status"] == "PASS" else 1)
