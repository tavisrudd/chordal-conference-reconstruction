#!/usr/bin/env python3
"""Exact Coxeter-number conic-phase certificate.

The symbolic path constructs the integral A3=D3 and B3 reflection
arrangements, their projective rank-two flats, and the secants of every
singular stratum.  The finite-field path independently enumerates PG(2,p)
for representative odd good primes.  H3's already-certified mirror and
source-orbit counts are read from the tracked good-reduction and arithmetic-phase certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
from typing import TypeAlias


SCHEMA = "coxeter-coxeter-number-conic-phase-v1"
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path(__file__).with_suffix(".json")
GOOD_REDUCTION = HERE / "h3_good_reduction.json"
ARITHMETIC_PHASE = HERE / "h3_arithmetic_phase.json"
DEEP_HOLE = HERE / "deep_hole_transform.json"
DUALITY = HERE / "common_duality.json"
DECORATED = HERE / "decorated_parent.json"

Vector: TypeAlias = tuple[int, int, int]
Matrix: TypeAlias = tuple[Vector, Vector, Vector]
Permutation: TypeAlias = tuple[int, ...]


def dot(left: Vector, right: Vector) -> int:
    return sum(a * b for a, b in zip(left, right))


def cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def det(rows: tuple[Vector, Vector, Vector]) -> int:
    return dot(rows[0], cross(rows[1], rows[2]))


def matrix_vector(matrix: Matrix, vector: Vector, prime: int) -> Vector:
    return tuple(dot(row, vector) % prime for row in matrix)  # type: ignore[return-value]


def matrix_product(left: Matrix, right: Matrix, prime: int) -> Matrix:
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(row[k] * column[k] for k in range(3)) % prime for column in columns)
        for row in left
    )  # type: ignore[return-value]


def normalize_matrix(matrix: Matrix, prime: int) -> Matrix:
    for row in matrix:
        for entry in row:
            if entry % prime:
                inverse = pow(entry % prime, -1, prime)
                return tuple(
                    tuple(value * inverse % prime for value in current) for current in matrix
                )  # type: ignore[return-value]
    raise AssertionError("zero matrix")


def compose_permutations(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse_permutation(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def generated_permutation_group(generators: list[Permutation]) -> set[Permutation]:
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            child = compose_permutations(current, generator)
            if child not in group:
                group.add(child)
                frontier.append(child)
    return group


def normalize_integer(vector: Vector) -> Vector:
    divisor = 0
    for entry in vector:
        divisor = gcd(divisor, abs(entry))
    assert divisor
    result = tuple(entry // divisor for entry in vector)
    for entry in result:
        if entry:
            return tuple(-value for value in result) if entry < 0 else result  # type: ignore[return-value]
    raise AssertionError("zero vector")


def normalize_mod(vector: Vector, prime: int) -> Vector:
    reduced = tuple(entry % prime for entry in vector)
    for entry in reduced:
        if entry:
            inverse = pow(entry, -1, prime)
            return tuple(value * inverse % prime for value in reduced)  # type: ignore[return-value]
    raise AssertionError("zero vector")


def arrangements() -> dict[str, list[Vector]]:
    long_roots = (
        [(1, sign, 0) for sign in (-1, 1)]
        + [(1, 0, sign) for sign in (-1, 1)]
        + [(0, 1, sign) for sign in (-1, 1)]
    )
    return {
        "A3": long_roots,
        "B3": [(1, 0, 0), (0, 1, 0), (0, 0, 1)] + long_roots,
    }


def singular_strata(roots: list[Vector]) -> dict[int, list[Vector]]:
    points = {
        normalize_integer(cross(left, right))
        for left, right in combinations(roots, 2)
        if cross(left, right) != (0, 0, 0)
    }
    strata: dict[int, list[Vector]] = {}
    for point in sorted(points):
        multiplicity = sum(dot(root, point) == 0 for root in roots)
        strata.setdefault(multiplicity, []).append(point)
    return strata


def secants(points: list[Vector]) -> set[Vector]:
    return {
        normalize_integer(cross(left, right))
        for left, right in combinations(points, 2)
        if cross(left, right) != (0, 0, 0)
    }


def stratum_certificate(roots: list[Vector], points: list[Vector]) -> dict[str, object]:
    triple_determinants = [det(triple) for triple in combinations(points, 3)]
    lines = secants(points)
    mirror_set = {normalize_integer(root) for root in roots}
    return {
        "point_count": len(points),
        "points": [list(point) for point in points],
        "arc": len(points) >= 3 and all(value != 0 for value in triple_determinants),
        "triple_determinant_absolute_counts": {
            str(key): value
            for key, value in sorted(Counter(abs(value) for value in triple_determinants).items())
        },
        "distinct_secant_count": len(lines),
        "secants_contained_in_mirrors": len(lines & mirror_set),
        "all_mirrors_are_secants": lines == mirror_set,
    }


def special_line_certificate(roots: list[Vector]) -> dict[str, object]:
    strata = singular_strata(roots)
    points = [point for values in strata.values() for point in values]
    multiplicity = {
        point: sum(dot(root, point) == 0 for root in roots) for point in points
    }
    mirrors = {normalize_integer(root) for root in roots}
    lines = {
        normalize_integer(cross(left, right))
        for left, right in combinations(points, 2)
        if cross(left, right) != (0, 0, 0)
    } - mirrors
    delta_counts: Counter[int] = Counter()
    special_through: dict[Vector, int] = {point: 0 for point in points}
    for line in lines:
        incident = [point for point in points if dot(line, point) == 0]
        delta_counts[sum(multiplicity[point] - 1 for point in incident)] += 1
        for point in incident:
            special_through[point] += 1
    pencil_types = Counter(
        (multiplicity[point], special_through[point]) for point in points
    )
    return {
        "special_nonmirror_line_count": len(lines),
        "special_nonmirror_delta_counts": {
            str(key): value for key, value in sorted(delta_counts.items())
        },
        "singular_point_pencil_types": [
            {
                "arrangement_multiplicity": key[0],
                "special_nonmirror_lines_through_point": key[1],
                "point_count": value,
            }
            for key, value in sorted(pencil_types.items())
        ],
    }


def triangular_source(mirror_count: int) -> dict[str, object]:
    discriminant = 1 + 8 * mirror_count
    root = isqrt(discriminant)
    source_size = (1 + root) // 2 if root * root == discriminant and root % 2 == 1 else None
    return {
        "mirror_count": mirror_count,
        "discriminant_1_plus_8m": discriminant,
        "arc_source_size_if_all_mirrors_are_secants": source_size,
    }


def symbolic_certificate() -> dict[str, object]:
    result: dict[str, object] = {}
    for name, roots in arrangements().items():
        strata = singular_strata(roots)
        nonzero_minors = [abs(det(triple)) for triple in combinations(roots, 3) if det(triple)]
        result[name] = {
            "mirror_count": len(roots),
            "rank_two_flat_multiplicity_counts": {
                str(multiplicity): len(points) for multiplicity, points in sorted(strata.items())
            },
            "nonzero_mirror_minor_absolute_counts": {
                str(key): value for key, value in sorted(Counter(nonzero_minors).items())
            },
            "bad_lattice_characteristics": [2],
            "singular_strata": {
                str(multiplicity): stratum_certificate(roots, points)
                for multiplicity, points in sorted(strata.items())
            },
            "special_line_ledger": special_line_certificate(roots),
        }

    a3 = result["A3"]
    b3 = result["B3"]
    assert isinstance(a3, dict) and isinstance(b3, dict)
    assert a3["rank_two_flat_multiplicity_counts"] == {"2": 3, "3": 4}
    assert b3["rank_two_flat_multiplicity_counts"] == {"2": 6, "3": 4, "4": 3}
    assert a3["bad_lattice_characteristics"] == [2]
    assert b3["bad_lattice_characteristics"] == [2]
    assert a3["singular_strata"]["3"]["all_mirrors_are_secants"]  # type: ignore[index]
    assert b3["singular_strata"]["3"]["arc"]  # type: ignore[index]
    assert b3["singular_strata"]["3"]["distinct_secant_count"] == 6  # type: ignore[index]
    assert b3["singular_strata"]["4"]["distinct_secant_count"] == 3  # type: ignore[index]
    assert a3["special_line_ledger"] == {  # type: ignore[comparison-overlap]
        "special_nonmirror_line_count": 3,
        "special_nonmirror_delta_counts": {"2": 3},
        "singular_point_pencil_types": [
            {
                "arrangement_multiplicity": 2,
                "special_nonmirror_lines_through_point": 2,
                "point_count": 3,
            },
            {
                "arrangement_multiplicity": 3,
                "special_nonmirror_lines_through_point": 0,
                "point_count": 4,
            },
        ],
    }
    assert b3["special_line_ledger"] == {  # type: ignore[comparison-overlap]
        "special_nonmirror_line_count": 16,
        "special_nonmirror_delta_counts": {"3": 16},
        "singular_point_pencil_types": [
            {
                "arrangement_multiplicity": 2,
                "special_nonmirror_lines_through_point": 4,
                "point_count": 6,
            },
            {
                "arrangement_multiplicity": 3,
                "special_nonmirror_lines_through_point": 3,
                "point_count": 4,
            },
            {
                "arrangement_multiplicity": 4,
                "special_nonmirror_lines_through_point": 0,
                "point_count": 3,
            },
        ],
    }

    good_reduction = json.loads(GOOD_REDUCTION.read_text())
    arithmetic_phase = json.loads(ARITHMETIC_PHASE.read_text())
    deep_hole = json.loads(DEEP_HOLE.read_text())
    h3_mirrors = sum(int(value) for value in good_reduction["symbolic"]["rank_2_flat_multiplicity_counts"].values())
    h3_source = good_reduction["symbolic"]["rank_2_flat_multiplicity_counts"]["5"]
    h3_secants = arithmetic_phase["characteristic_5_source_phase"]["secant_count"]
    assert (h3_mirrors, h3_source, h3_secants) == (31, 6, 15)
    assert arithmetic_phase["q11_deep_hole_phase"]["field_order"] == 11
    assert len(arithmetic_phase["q11_deep_hole_phase"]["deep_hole_points_h3_coordinates"]) == 12
    assert arithmetic_phase["q11_deep_hole_phase"]["child_projective_code_parameters"] == [12, 3, 10]
    # There are 31 singular flats, but the reflection arrangement has 15 mirrors.
    h3_reflection_mirrors = good_reduction["finite_field_replays"][1]["distinct_mirrors"]
    assert h3_reflection_mirrors == 15

    phase_table = {
        "A3": {
            "coxeter_number_h": 4,
            "middle_exponent_e": 2,
            "mirror_count_3h_over_2": 6,
            "complement_length": "(q-2)(q-3)",
            "line_delta_counts": {"0": "(q-3)^2", "1": "3(q-3)", "2": "4q-5"},
            "line_intersection_sizes": {"0": "q-5", "1": "q-4", "2": "q-3"},
            "maximum_line_intersection": "q-3",
            "evaluation_code_distance": "(q-3)^2",
            "coxeter_conic_field": 5,
        },
        "B3": {
            "coxeter_number_h": 6,
            "middle_exponent_e": 3,
            "mirror_count_3h_over_2": 9,
            "complement_length": "(q-3)(q-5)",
            "line_delta_counts": {
                "0": "(q-5)(q-7)",
                "1": "6(q-5)",
                "2": "4(q-5)",
                "3": "3q+7",
            },
            "line_intersection_sizes": {
                "0": "q-8",
                "1": "q-7",
                "2": "q-6",
                "3": "q-5",
            },
            "maximum_line_intersection": "q-5",
            "evaluation_code_distance": "(q-4)(q-5)",
            "coxeter_conic_field": 7,
        },
        "H3": {
            "coxeter_number_h": 10,
            "middle_exponent_e": 5,
            "mirror_count_3h_over_2": 15,
            "complement_length": "(q-5)(q-9)",
            "line_delta_counts": deep_hole["all_nonmirror_line_delta_counts"],
            "line_intersection_sizes": deep_hole["intersection_sizes"],
            "maximum_line_intersection": "q-9",
            "evaluation_code_distance": "(q-6)(q-9)",
            "coxeter_conic_field": 11,
        },
    }
    for data in phase_table.values():
        h = data["coxeter_number_h"]
        e = data["middle_exponent_e"]
        assert h == 2 * e
        assert data["mirror_count_3h_over_2"] == 3 * h // 2
        assert data["coxeter_conic_field"] == h + 1

    obstruction = {
        "A3": triangular_source(6),
        "B3": triangular_source(9),
        "H3": triangular_source(15),
    }
    assert obstruction["A3"]["arc_source_size_if_all_mirrors_are_secants"] == 4
    assert obstruction["B3"]["arc_source_size_if_all_mirrors_are_secants"] is None
    assert obstruction["H3"]["arc_source_size_if_all_mirrors_are_secants"] == 6
    return {
        "integral_arrangements": result,
        "h3_dependency": {
            "reflection_mirror_count": h3_reflection_mirrors,
            "fivefold_source_orbit_size": h3_source,
            "source_secant_count": h3_secants,
            "good_reduction_sha256": hashlib.sha256(GOOD_REDUCTION.read_bytes()).hexdigest(),
            "arithmetic_phase_sha256": hashlib.sha256(ARITHMETIC_PHASE.read_bytes()).hexdigest(),
            "deep_hole_sha256": hashlib.sha256(DEEP_HOLE.read_bytes()).hexdigest(),
        },
        "coxeter_number_phase_table": phase_table,
        "common_theorem": {
            "length": "n_T(q)=(q-h/2)(q-h+1)",
            "distance": "d_T(q)=(q-h/2-1)(q-h+1)",
            "coxeter_conic_phase": "q=h+1 gives the full invariant conic and [q+1,3,q-1] GRS",
            "stable_recovery": "q>3h/2-1 makes the reflection mirrors exactly the disjoint lines",
        },
        "secant_count_obstruction": obstruction,
        "b3_defect_verdict": (
            "No full A3/B3/H3 complete-deep-hole functor exists with an MDS "
            "source arc whose complete secants are exactly the reflection mirrors: B3 has nine "
            "mirrors, while an n-arc has exactly binomial(n,2) distinct secants. Its natural "
            "four-point singular orbit sees only the six-mirror long-root D3=A3 subsystem, and "
            "its three-point orbit sees only the three coordinate mirrors."
        ),
    }


def projective_points(prime: int) -> list[Vector]:
    points: list[Vector] = []
    for x in range(prime):
        for y in range(prime):
            for z in range(prime):
                if (x, y, z) != (0, 0, 0):
                    point = normalize_mod((x, y, z), prime)
                    if point not in points:
                        points.append(point)
    return sorted(points)


def normalize_pair(pair: tuple[int, int], prime: int) -> tuple[int, int]:
    for entry in pair:
        if entry % prime:
            inverse = pow(entry % prime, -1, prime)
            return pair[0] * inverse % prime, pair[1] * inverse % prime
    raise AssertionError("zero pair")


def conic_parameterization(prime: int) -> tuple[list[Vector], list[tuple[int, int]]]:
    conic = [point for point in projective_points(prime) if dot(point, point) % prime == 0]
    base = conic[0]
    pencil = [line for line in projective_points(prime) if dot(line, base) % prime == 0]
    first = pencil[0]
    second = next(line for line in pencil[1:] if cross(first, line) != (0, 0, 0))
    parameters = [normalize_pair(pair, prime) for pair in (
        [(1, value) for value in range(prime)] + [(0, 1)]
    )]
    points_by_parameter: list[Vector] = []
    for left, right in parameters:
        line = normalize_mod(
            tuple((left * first[index] + right * second[index]) % prime for index in range(3)),
            prime,
        )
        incident = [point for point in conic if dot(line, point) % prime == 0]
        assert len(incident) in (1, 2) and base in incident
        points_by_parameter.append(base if len(incident) == 1 else next(p for p in incident if p != base))
    assert len(set(points_by_parameter)) == prime + 1
    return points_by_parameter, parameters


def reflection_matrix(root: Vector, prime: int) -> Matrix:
    length = dot(root, root) % prime
    matrix = tuple(
        tuple(
            (length * (1 if i == j else 0) - 2 * root[i] * root[j]) % prime
            for j in range(3)
        )
        for i in range(3)
    )
    assert det(matrix) % prime
    return normalize_matrix(matrix, prime)  # type: ignore[arg-type]


def conic_symmetry_certificate(name: str, prime: int, roots: list[Vector]) -> dict[str, object]:
    conic, parameters = conic_parameterization(prime)
    point_index = {point: index for index, point in enumerate(conic)}
    reflection_permutations: list[Permutation] = []
    for root in roots:
        matrix = reflection_matrix(root, prime)
        permutation = tuple(
            point_index[normalize_mod(matrix_vector(matrix, point, prime), prime)] for point in conic
        )
        reflection_permutations.append(permutation)
    coxeter_group = generated_permutation_group(reflection_permutations)

    parameter_index = {parameter: index for index, parameter in enumerate(parameters)}
    normalized_2x2: set[tuple[int, int, int, int]] = set()
    for a in range(prime):
        for b in range(prime):
            for c in range(prime):
                for d in range(prime):
                    if (a * d - b * c) % prime == 0:
                        continue
                    entries = (a, b, c, d)
                    pivot = next(value for value in entries if value)
                    inverse = pow(pivot, -1, prime)
                    normalized_2x2.add(tuple(value * inverse % prime for value in entries))
    full_group: set[Permutation] = set()
    for a, b, c, d in normalized_2x2:
        full_group.add(tuple(
            parameter_index[
                normalize_pair((a * left + b * right, c * left + d * right), prime)
            ]
            for left, right in parameters
        ))
    assert len(full_group) == prime * (prime * prime - 1)
    assert coxeter_group <= full_group
    normalizer = {
        permutation
        for permutation in full_group
        if {
            compose_permutations(
                compose_permutations(permutation, element),
                inverse_permutation(permutation),
            )
            for element in coxeter_group
        } == coxeter_group
    }
    assert len({permutation[0] for permutation in coxeter_group}) == prime + 1
    expected_coxeter_order = 24
    assert len(coxeter_group) == expected_coxeter_order
    assert len(normalizer) == expected_coxeter_order
    return {
        "type": name,
        "field_order": prime,
        "conic_point_count": len(conic),
        "projective_coxeter_group_order": len(coxeter_group),
        "full_conic_stabilizer_order": len(full_group),
        "coxeter_group_transitive_on_conic": True,
        "coxeter_normalizer_order": len(normalizer),
        "conjugate_coxeter_decoration_count": len(full_group) // len(normalizer),
    }


def finite_replay(prime: int, name: str, roots: list[Vector]) -> dict[str, object]:
    mirrors = {normalize_mod(root, prime) for root in roots}
    points = projective_points(prime)
    spectrum = Counter(
        sum(dot(root, point) % prime == 0 for root in mirrors) for point in points
    )
    complement_points = {
        point
        for point in points
        if all(dot(root, point) % prime != 0 for root in mirrors)
    }
    complement = len(complement_points)
    expected = (prime - 2) * (prime - 3) if name == "A3" else (prime - 3) * (prime - 5)
    assert complement == expected
    line_intersections = Counter(
        sum(dot(line, point) % prime == 0 for point in complement_points)
        for line in points
    )
    maximum_line_intersection = max(line_intersections)
    spans_plane = any(det(triple) % prime != 0 for triple in combinations(complement_points, 3))
    if complement:
        expected_maximum = prime - 3 if name == "A3" else prime - 5
        assert maximum_line_intersection == expected_maximum
        assert spans_plane
        expected_distance = (
            (prime - 3) ** 2 if name == "A3" else (prime - 4) * (prime - 5)
        )
        assert complement - maximum_line_intersection == expected_distance
    conic_points = {point for point in points if dot(point, point) % prime == 0}
    coxeter_conic_field = 5 if name == "A3" else 7
    if prime == coxeter_conic_field:
        assert complement_points == conic_points
    return {
        "type": name,
        "characteristic": prime,
        "distinct_mirrors": len(mirrors),
        "projective_point_multiplicity_spectrum": {
            str(key): value for key, value in sorted(spectrum.items())
        },
        "arrangement_complement_size": complement,
        "expected_complement_polynomial_value": expected,
        "invariant_conic_size": len(conic_points),
        "complement_equals_invariant_conic": complement_points == conic_points,
        "line_intersection_spectrum": {
            str(key): value for key, value in sorted(line_intersections.items())
        },
        "maximum_line_intersection": maximum_line_intersection,
        "spans_projective_plane": spans_plane,
        "evaluation_code_parameters": (
            [complement, 3, complement - maximum_line_intersection] if complement else None
        ),
    }


def build_certificate() -> dict[str, object]:
    models = arrangements()
    duality = json.loads(DUALITY.read_text())
    decorated = json.loads(DECORATED.read_text())
    h3_symmetry = {
        "type": "H3",
        "field_order": 11,
        "conic_point_count": duality["conic_projective_point_count"],
        "projective_coxeter_group_order": decorated["marked_fibre"]["fixed_a5_order"],
        "full_conic_stabilizer_order": decorated["marked_fibre"]["full_child_stabilizer_order"],
        "coxeter_group_transitive_on_conic": True,
        "coxeter_normalizer_order": decorated["marked_fibre"]["fixed_a5_normalizer_order"],
        "conjugate_coxeter_decoration_count": decorated["marked_fibre"]["conjugate_a5_parent_count"],
        "duality_sha256": hashlib.sha256(DUALITY.read_bytes()).hexdigest(),
        "decorated_sha256": hashlib.sha256(DECORATED.read_bytes()).hexdigest(),
    }
    assert h3_symmetry["conic_point_count"] == 12
    assert h3_symmetry["full_conic_stabilizer_order"] == 1320
    assert h3_symmetry["conjugate_coxeter_decoration_count"] == 22
    rank_four_gate = []
    for name, h, q, exponents in (
        ("A4", 5, 6, (1, 2, 3, 4)),
        ("B4", 8, 9, (1, 3, 5, 7)),
        ("D4", 6, 7, (1, 3, 3, 5)),
        ("F4", 12, 13, (1, 5, 7, 11)),
        ("H4", 30, 31, (1, 11, 19, 29)),
    ):
        field_exists = q in (7, 9, 13, 31)
        complement = None
        quadric_sizes = None
        if field_exists:
            complement = 1
            for exponent in exponents[1:]:
                complement *= q - exponent
            quadric_sizes = sorted((q * q + 1, (q + 1) * (q + 1)))
            assert complement not in quadric_sizes
        rank_four_gate.append({
            "type": name,
            "coxeter_number_h": h,
            "candidate_q_h_plus_1": q,
            "q_is_prime_power": field_exists,
            "projective_complement_size": complement,
            "possible_nonsingular_quadric_surface_sizes": quadric_sizes,
            "quadric_equality": False,
        })
    return {
        "schema": SCHEMA,
        "symbolic": symbolic_certificate(),
        "independent_finite_field_replays": [
            finite_replay(prime, name, models[name])
            for prime in (3, 5, 7, 11)
            for name in ("A3", "B3")
        ],
        "free_symmetry_completion_upgrade": [
            conic_symmetry_certificate("A3", 5, models["A3"]),
            conic_symmetry_certificate("B3", 7, models["B3"]),
            h3_symmetry,
        ],
        "first_higher_rank_gate": {
            "rank": 4,
            "types": rank_four_gate,
            "verdict": "No irreducible rank-four Coxeter type has the q=h+1 complement-count equality with a nonsingular quadric surface.",
        },
    }


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("choose exactly one of --write or --check")
    expected = canonical_bytes(build_certificate())
    if args.write:
        DEFAULT_OUTPUT.write_bytes(expected)
        print(f"wrote {DEFAULT_OUTPUT.name}: {len(expected)} bytes")
        return
    actual = DEFAULT_OUTPUT.read_bytes()
    assert actual == expected, f"stale certificate: regenerate {DEFAULT_OUTPUT.name}"
    print(
        f"ok {DEFAULT_OUTPUT.name}: {len(actual)} bytes, "
        f"sha256={hashlib.sha256(actual).hexdigest()}"
    )


if __name__ == "__main__":
    main()
