#!/usr/bin/env python3
"""Exact reconnaissance for the Paper-II-to-six-axis cubic bridge.

This script uses Paper II's frozen matching arithmetic but computes the new
restriction independently: the H3 quotient is projected to its unique
five-dimensional A5 constituent and compared with every A5-invariant
order-six conference triangle cubic over F_11.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
FROZEN = HERE / "frozen"
MATCHING_PATH = FROZEN / "matching_module.py"
SCOUT_PATH = FROZEN / "matching_orbit_scout.json"
INPUT_PATHS = tuple(
    FROZEN / name
    for name in (
        "matching_module.py",
        "matching_orbit_scout.json",
        "coxeter_conic_phase.py",
        "common_duality.py",
        "common_duality.json",
        "a5_subgroup_decoder.py",
        "h3_good_reduction.json",
        "h3_arithmetic_phase.json",
        "deep_hole_transform.json",
        "decorated_parent.json",
        "decorated_parent_replay.py",
    )
)
OUTPUT = HERE / "paper_ii_chordal_axis.json"
SCHEMA = "paper-v-chordal-axis-v2"
P = 11


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("paper_v_matching", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_module(MATCHING_PATH)


def columns(matrix):
    return [list(column) for column in zip(*matrix)]


def mat_add(left, right):
    return [
        [(a + b) % P for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def mat_scale(scalar, matrix):
    return [[scalar * value % P for value in row] for row in matrix]


def projective_normalize(vector):
    pivot = next(value for value in vector if value % P)
    inverse = pow(pivot, -1, P)
    return [inverse * value % P for value in vector]


def matrix_rank(matrix):
    return M.rank(matrix, P)


def h3_data():
    scout = json.loads(SCOUT_PATH.read_text())
    record = next(item for item in scout["types"] if item["type"] == "H3")
    conic, parameters = M.COXETER.conic_parameterization(P)
    endpoints = tuple(parameters)
    full_group, psl_group = M.full_pgl(P, parameters)
    parent_group = M.h3_group(P, conic)
    base_matching = tuple(tuple(pair) for pair in record["coxeter_invariant_matching"])
    assert all(M.matching_image(element, base_matching) == base_matching for element in parent_group)
    orbit = sorted({M.matching_image(element, base_matching) for element in full_group})
    orbit_index = {matching: index for index, matching in enumerate(orbit)}

    base_product = M.matching_product(base_matching, endpoints, P)
    quotient_vectors = []
    for matching in orbit:
        product = M.matching_product(matching, endpoints, P)
        difference = {
            exponent: (product.get(exponent, 0) - base_product.get(exponent, 0)) % P
            for exponent in set(product) | set(base_product)
        }
        quotient_vectors.append(M.quotient_by_conic(difference, 4, P))

    image_matrix = M.transpose(quotient_vectors)
    _reduced, coordinate_pivots = M.rref(M.transpose(image_matrix), P)
    assert len(coordinate_pivots) == 10
    reduced_vectors = [
        [vector[index] for index in coordinate_pivots] for vector in quotient_vectors
    ]
    _point_reduced, point_basis_indices = M.rref(M.transpose(reduced_vectors), P)
    assert len(point_basis_indices) == 10
    point_basis = M.transpose([reduced_vectors[index] for index in point_basis_indices])
    point_basis_inverse = M.matrix_inverse(point_basis, P)

    def induced_action(element):
        action = M.action_permutation(element, orbit, orbit_index)
        target_basis = M.transpose([reduced_vectors[action[index]] for index in point_basis_indices])
        result = M.matrix_product(target_basis, point_basis_inverse, P)
        assert all(
            M.matrix_vector(result, reduced_vectors[index], P) == reduced_vectors[action[index]]
            for index in range(len(orbit))
        )
        return result

    classes = M.labelled_classes("H3", parent_group)
    class_by_element = {
        element: label for label, _representative, elements in classes for element in elements
    }
    actions = {element: induced_action(element) for element in parent_group}
    character = dict(zip(M.character_table("H3", P)[0], M.character_table("H3", P)[1]["5"]))
    projector = [[0] * 10 for _ in range(10)]
    scalar = 5 * pow(len(parent_group), -1, P) % P
    for element in parent_group:
        projector = mat_add(
            projector,
            mat_scale(scalar * character[class_by_element[element]] % P, actions[element]),
        )
    assert M.matrix_product(projector, projector, P) == projector
    assert matrix_rank(projector) == 5

    pivot_columns = M.rref(projector, P)[1]
    assert len(pivot_columns) == 5
    constituent_basis = [[projector[row][column] for column in pivot_columns] for row in range(10)]
    pivot_rows = M.rref(M.transpose(constituent_basis), P)[1]
    assert len(pivot_rows) == 5
    coordinate_inverse = M.matrix_inverse(
        [[constituent_basis[row][column] for column in range(5)] for row in pivot_rows], P
    )

    def constituent_coordinates(vector):
        projected = M.matrix_vector(projector, vector, P)
        return M.matrix_vector(coordinate_inverse, [projected[row] for row in pivot_rows], P)

    constituent_actions = {}
    for element in parent_group:
        moved_basis = M.matrix_product(actions[element], constituent_basis, P)
        constituent_actions[element] = M.matrix_product(
            coordinate_inverse,
            [[moved_basis[row][column] for column in range(5)] for row in pivot_rows],
            P,
        )

    base_pairs = list(base_matching)
    pair_index = {pair: index for index, pair in enumerate(base_pairs)}
    pair_permutations = {}
    axis_actions = {}
    for element in parent_group:
        permutation = tuple(
            pair_index[tuple(sorted((element[left], element[right])))]
            for left, right in base_pairs
        )
        pair_permutations[element] = permutation
        matrix = [[0] * 5 for _ in range(5)]
        for source in range(5):
            image_vector = [0] * 6
            image_vector[permutation[source]] = 1
            image_vector[permutation[5]] = -1 % P
            for row in range(5):
                matrix[row][source] = image_vector[row]
        axis_actions[element] = matrix

    generators = M.permutation_generators(parent_group)
    equations = []
    for element in generators:
        left_action = constituent_actions[element]
        right_action = axis_actions[element]
        for row in range(5):
            for column in range(5):
                equation = [0] * 25
                for middle in range(5):
                    equation[middle * 5 + column] = (
                        equation[middle * 5 + column] + left_action[row][middle]
                    ) % P
                    equation[row * 5 + middle] = (
                        equation[row * 5 + middle] - right_action[middle][column]
                    ) % P
                equations.append(equation)
    intertwiners = M.nullspace(equations, P)
    assert len(intertwiners) == 1
    intertwiner = [intertwiners[0][row * 5 : (row + 1) * 5] for row in range(5)]
    assert matrix_rank(intertwiner) == 5
    inverse_intertwiner = M.matrix_inverse(intertwiner, P)

    def axis_coordinates(vector):
        return M.matrix_vector(inverse_intertwiner, constituent_coordinates(vector), P)

    sheets = []
    unseen = set(orbit)
    while unseen:
        representative = min(unseen)
        sheet = {M.matching_image(element, representative) for element in psl_group}
        unseen -= sheet
        sheets.append(sheet)
    sheets.sort(key=lambda sheet: sorted(sheet))
    assert len(sheets) == 2
    signs = [1 if matching in sheets[0] else -1 % P for matching in orbit]

    gram = [[(1 if row == column else 0) + 1 for column in range(5)] for row in range(5)]
    covectors = [M.matrix_vector(gram, axis_coordinates(vector), P) for vector in reduced_vectors]
    tensor = [
        sum(sign * monomial for sign, monomial in zip(signs, values)) % P
        for values in zip(*(M.symmetric_power(vector, 3, P) for vector in covectors))
    ]
    cube_basis = list(itertools.combinations_with_replacement(range(5), 3))
    cubic_polynomial = []
    for coefficient, indices in zip(tensor, cube_basis):
        multiplicity = 6
        if indices[0] == indices[2]:
            multiplicity = 1
        elif indices[0] == indices[1] or indices[1] == indices[2]:
            multiplicity = 3
        cubic_polynomial.append(multiplicity * coefficient % P)

    normalization_pivot = next(value for value in cubic_polynomial if value % P)

    normalized_cubic = projective_normalize(cubic_polynomial)
    assert cubic_polynomial == [
        normalization_pivot * value % P for value in normalized_cubic
    ]

    return {
        "parent_group": parent_group,
        "pair_permutations": pair_permutations,
        "axis_actions": axis_actions,
        "parameters": endpoints,
        "base_matching": base_matching,
        "intertwiner": intertwiner,
        "generator_permutations": [list(element) for element in generators],
        "generator_axis_actions": [axis_actions[element] for element in generators],
        "generator_constituent_actions": [
            constituent_actions[element] for element in generators
        ],
        "projected_cubic_raw": cubic_polynomial,
        "projected_cubic_normalization_pivot": normalization_pivot,
        "projected_cubic": normalized_cubic,
        "hom_dimension": len(intertwiners),
        "projector_rank": matrix_rank(projector),
    }


def conference_candidates(data):
    edges = list(itertools.combinations(range(1, 6), 2))
    basis = list(itertools.combinations_with_replacement(range(5), 3))
    basis_index = {indices: index for index, indices in enumerate(basis)}
    augmentation_rows = [
        [((1 if axis == coordinate else 0) - (1 if axis == 5 else 0)) % P for coordinate in range(5)]
        for axis in range(6)
    ]

    def polynomial(matrix):
        result = [0] * len(basis)
        for i, j, k in itertools.combinations(range(6), 3):
            tau = matrix[i][j] * matrix[j][k] * matrix[k][i]
            for a in range(5):
                for b in range(5):
                    for c in range(5):
                        index = basis_index[tuple(sorted((a, b, c)))]
                        result[index] = (
                            result[index]
                            + tau
                            * augmentation_rows[i][a]
                            * augmentation_rows[j][b]
                            * augmentation_rows[k][c]
                        ) % P
        return projective_normalize(result)

    records = []
    for bits in itertools.product((-1, 1), repeat=10):
        matrix = [[0] * 6 for _ in range(6)]
        for index in range(1, 6):
            matrix[0][index] = matrix[index][0] = 1
        for (left, right), sign in zip(edges, bits):
            matrix[left][right] = matrix[right][left] = sign
        square = M.matrix_product(matrix, matrix, P)
        if square != [[5 if row == column else 0 for column in range(6)] for row in range(6)]:
            continue
        triangle = {
            triple: matrix[triple[0]][triple[1]]
            * matrix[triple[1]][triple[2]]
            * matrix[triple[2]][triple[0]]
            for triple in itertools.combinations(range(6), 3)
        }
        invariant = all(
            all(
                triangle[tuple(sorted(permutation[index] for index in triple))] == triangle[triple]
                for triple in triangle
            )
            for permutation in data["pair_permutations"].values()
        )
        records.append({"matrix": matrix, "invariant": invariant, "polynomial": polynomial(matrix)})
    assert len(records) == 12
    return records


def file_record(path):
    payload = path.read_bytes()
    return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def projective_singular_points(coefficients):
    basis = list(itertools.combinations_with_replacement(range(5), 3))
    singular = []
    for pivot in range(5):
        for tail in itertools.product(range(P), repeat=4 - pivot):
            point = (0,) * pivot + (1,) + tail
            gradient = [0] * 5
            for coefficient, indices in zip(coefficients, basis):
                for position, index in enumerate(indices):
                    product = coefficient
                    for other_position, other_index in enumerate(indices):
                        if other_position != position:
                            product = product * point[other_index] % P
                    gradient[index] = (gradient[index] + product) % P
            if not any(gradient):
                singular.append(point)
    return singular


def normalize_point(point):
    return tuple(projective_normalize(list(point)))


def orbit_record(points, data):
    point_set = set(points)
    stabilizers = []
    for point in points:
        stabilizer = frozenset(
            element
            for element in data["parent_group"]
            if normalize_point(M.matrix_vector(data["axis_actions"][element], point, P)) == point
        )
        stabilizers.append(stabilizer)
    assert all(
        normalize_point(M.matrix_vector(action, point, P)) in point_set
        for action in data["axis_actions"].values()
        for point in points
    )
    fibres = {}
    for index, stabilizer in enumerate(stabilizers):
        fibres.setdefault(stabilizer, []).append(index)
    return {
        "point_count": len(points),
        "stabilizer_orders": sorted({len(stabilizer) for stabilizer in stabilizers}),
        "equal_stabilizer_fibre_sizes": sorted(len(fibre) for fibre in fibres.values()),
        "transitive": len(
            {
                normalize_point(M.matrix_vector(action, points[0], P))
                for action in data["axis_actions"].values()
            }
        )
        == len(points),
    }


def transform_cubic(coefficients, matrix):
    basis = list(itertools.combinations_with_replacement(range(5), 3))
    basis_index = {indices: index for index, indices in enumerate(basis)}
    result = [0] * len(basis)
    for coefficient, (i, j, k) in zip(coefficients, basis):
        for left in range(5):
            for middle in range(5):
                for right in range(5):
                    index = basis_index[tuple(sorted((left, middle, right)))]
                    result[index] = (
                        result[index]
                        + coefficient
                        * matrix[i][left]
                        * matrix[j][middle]
                        * matrix[k][right]
                    ) % P
    return result


def outer_pencil_action(data, conference, chordal):
    """Return the exact outer-normalizer action on the invariant pencil."""
    pair_group = set(data["pair_permutations"].values())

    def compose(left, right):
        return tuple(left[right[index]] for index in range(6))

    def invert(permutation):
        result = [0] * 6
        for source, image in enumerate(permutation):
            result[image] = source
        return tuple(result)

    normalizer = []
    for permutation in itertools.permutations(range(6)):
        inverse = invert(permutation)
        if {
            compose(compose(permutation, element), inverse) for element in pair_group
        } == pair_group:
            normalizer.append(permutation)
    assert len(normalizer) == 120
    outer = next(permutation for permutation in normalizer if permutation not in pair_group)

    action = [[0] * 5 for _ in range(5)]
    for source in range(5):
        image_vector = [0] * 6
        image_vector[outer[source]] = 1
        image_vector[outer[5]] = -1 % P
        for row in range(5):
            action[row][source] = image_vector[row]

    def pencil_coordinates(cubic):
        for left in range(P):
            for right in range(P):
                if not (left or right):
                    continue
                if all(
                    (left * c + right * h - value) % P == 0
                    for c, h, value in zip(conference, chordal, cubic)
                ):
                    return [left, right]
        raise AssertionError("outer image left the invariant pencil")

    conference_image = pencil_coordinates(transform_cubic(conference, action))
    chordal_image = pencil_coordinates(transform_cubic(chordal, action))
    assert conference_image == [P - 1, 0]
    assert chordal_image == [8, 1]
    return {
        "outer_axis_permutation": list(outer),
        "basis": ["conference", "paper_ii_chordal"],
        "conference_image_coordinates": conference_image,
        "paper_ii_chordal_image_coordinates": chordal_image,
        "projective_parameter_action": "[a:b] -> [-a+8b:b]",
        "fixed_conference_parameter": [1, 0],
        "fixed_ten_point_parameter": [1, 3],
        "exchanged_chordal_parameters": [[0, 1], [1, 7]],
    }


def chordal_certificate(data, cubic):
    singular = projective_singular_points(cubic)
    assert len(singular) == 12
    pair_group = set(data["pair_permutations"].values())
    element_by_pair_action = {
        action: element for element, action in data["pair_permutations"].items()
    }

    def compose(left, right):
        return tuple(left[right[index]] for index in range(6))

    def invert(permutation):
        result = [0] * 6
        for source, image in enumerate(permutation):
            result[image] = source
        return tuple(result)

    normalizer = []
    for permutation in itertools.permutations(range(6)):
        inverse = invert(permutation)
        if {
            compose(compose(permutation, element), inverse) for element in pair_group
        } == pair_group:
            normalizer.append(permutation)
    assert len(normalizer) == 120
    outer = next(permutation for permutation in normalizer if permutation not in pair_group)
    outer_inverse = invert(outer)
    outer_automorphism = {
        element: element_by_pair_action[
            compose(
                compose(outer, data["pair_permutations"][element]), outer_inverse
            )
        ]
        for element in data["parent_group"]
    }

    def mobius_matrix(element):
        equations = []
        for endpoint, (s, t) in enumerate(data["parameters"]):
            image_s, image_t = data["parameters"][element[endpoint]]
            first = [0] * 16
            first[0], first[1], first[4 + endpoint] = s, t, -image_s % P
            second = [0] * 16
            second[2], second[3], second[4 + endpoint] = s, t, -image_t % P
            equations.extend((first, second))
        solutions = M.nullspace(equations, P)
        assert len(solutions) == 1
        a, b, c, d = solutions[0][:4]
        determinant = (a * d - b * c) % P
        determinant_scale = next(
            value for value in range(1, P) if value * value * determinant % P == 1
        )
        return [
            [a * determinant_scale % P, b * determinant_scale % P],
            [c * determinant_scale % P, d * determinant_scale % P],
        ]

    def polynomial_product(left, right):
        result = [0] * (len(left) + len(right) - 1)
        for left_index, left_value in enumerate(left):
            for right_index, right_value in enumerate(right):
                result[left_index + right_index] = (
                    result[left_index + right_index] + left_value * right_value
                ) % P
        return result

    def polynomial_power(polynomial, exponent):
        result = [1]
        for _ in range(exponent):
            result = polynomial_product(result, polynomial)
        return result

    def fourth_symmetric_power(matrix):
        return [
            polynomial_product(
                polynomial_power(matrix[0], 4 - exponent),
                polynomial_power(matrix[1], exponent),
            )
            for exponent in range(5)
        ]

    def symmetric_fourth_intertwiners(twist):
        equations = []
        for element in M.permutation_generators(data["parent_group"]):
            left_element = outer_automorphism[element] if twist else element
            left_action = data["axis_actions"][left_element]
            right_action = fourth_symmetric_power(mobius_matrix(element))
            for row in range(5):
                for column in range(5):
                    equation = [0] * 25
                    for middle in range(5):
                        equation[middle * 5 + column] = (
                            equation[middle * 5 + column] + left_action[row][middle]
                        ) % P
                        equation[row * 5 + middle] = (
                            equation[row * 5 + middle] - right_action[middle][column]
                        ) % P
                    equations.append(equation)
        return M.nullspace(equations, P)

    direct_intertwiners = symmetric_fourth_intertwiners(False)
    intertwiners = symmetric_fourth_intertwiners(True)
    assert len(direct_intertwiners) == 1
    assert len(intertwiners) == 1
    direct_projectivity = [
        direct_intertwiners[0][row * 5 : (row + 1) * 5] for row in range(5)
    ]
    projectivity = [
        intertwiners[0][row * 5 : (row + 1) * 5] for row in range(5)
    ]
    assert matrix_rank(direct_projectivity) == 5
    assert matrix_rank(projectivity) == 5

    curve_points = []
    for s, t in data["parameters"]:
        veronese = [pow(s, 4 - exponent, P) * pow(t, exponent, P) % P for exponent in range(5)]
        curve_points.append(
            normalize_point(M.matrix_vector(projectivity, veronese, P))
        )
    assert set(curve_points) == set(singular)

    direct_curve_points = []
    for s, t in data["parameters"]:
        veronese = [pow(s, 4 - exponent, P) * pow(t, exponent, P) % P for exponent in range(5)]
        direct_curve_points.append(
            normalize_point(M.matrix_vector(direct_projectivity, veronese, P))
        )
    direct_curve_match = set(direct_curve_points) == set(singular)
    assert not direct_curve_match

    basis = list(itertools.combinations_with_replacement(range(5), 3))
    basis_index = {indices: index for index, indices in enumerate(basis)}
    hankel = [0] * len(basis)
    for indices, coefficient in {
        (0, 2, 4): 1,
        (1, 2, 3): 2,
        (0, 3, 3): -1,
        (1, 1, 4): -1,
        (2, 2, 2): -1,
    }.items():
        hankel[basis_index[indices]] = coefficient % P

    transformed = transform_cubic(cubic, projectivity)
    direct_transformed = transform_cubic(cubic, direct_projectivity)
    assert projective_normalize(transformed) == projective_normalize(hankel)
    direct_hankel_match = projective_normalize(direct_transformed) == projective_normalize(hankel)
    assert not direct_hankel_match
    scale = next(
        transformed[index] * pow(hankel[index], -1, P) % P
        for index in range(len(hankel))
        if hankel[index]
    )
    assert transformed == [scale * coefficient % P for coefficient in hankel]
    return {
        "singular_curve": "rational normal quartic",
        "f11_rational_points": 12,
        "axis_action_normalizer_order": len(normalizer),
        "direct_intertwiner_dimension": len(direct_intertwiners),
        "direct_veronese_curve_match": direct_curve_match,
        "direct_hankel_match": direct_hankel_match,
        "outer_twisted_intertwiner_dimension": len(intertwiners),
        "projectivity": projectivity,
        "hankel_scale": scale,
        "standard_equation": "det[[z0,z1,z2],[z1,z2,z3],[z2,z3,z4]]",
    }


def build_certificate():
    data = h3_data()
    candidates = conference_candidates(data)
    invariant = [record for record in candidates if record["invariant"]]
    matches = [record for record in invariant if record["polynomial"] == data["projected_cubic"]]
    distinct_lines = {tuple(record["polynomial"]) for record in invariant}
    conference_line = invariant[0]["polynomial"] if invariant else []
    projected_singular = projective_singular_points(data["projected_cubic"])
    conference_singular = projective_singular_points(conference_line)
    pencil = []
    for left, right in [(1, value) for value in range(P)] + [(0, 1)]:
        cubic = [
            (left * conference + right * projected) % P
            for conference, projected in zip(conference_line, data["projected_cubic"])
        ]
        singular = projective_singular_points(cubic)
        pencil.append(
            {
                "projective_parameter": [left, right],
                "f11_singular_point_count": len(singular),
                "f11_orbit_record": orbit_record(singular, data) if singular else None,
            }
        )
    return {
        "schema": SCHEMA,
        "verdict": "MATCH" if matches else "NO_MATCH",
        "field": 11,
        "base_matching": [list(pair) for pair in data["base_matching"]],
        "five_constituent_projector_rank": data["projector_rank"],
        "axis_to_quotient_hom_dimension": data["hom_dimension"],
        "intertwiner": data["intertwiner"],
        "generator_permutations": data["generator_permutations"],
        "generator_axis_actions": data["generator_axis_actions"],
        "generator_constituent_actions": data["generator_constituent_actions"],
        "projected_sheet_cubic": data["projected_cubic"],
        "projected_sheet_cubic_raw": data["projected_cubic_raw"],
        "projected_sheet_cubic_normalization_pivot": data[
            "projected_cubic_normalization_pivot"
        ],
        "projected_sheet_cubic_normalized_byte_sha256": hashlib.sha256(
            bytes(data["projected_cubic"])
        ).hexdigest(),
        "conference_triangle_cubic": conference_line,
        "projected_and_conference_span_dimension": matrix_rank(
            [data["projected_cubic"], conference_line]
        ),
        "conference_gauge_count": len(candidates),
        "a5_invariant_conference_gauge_count": len(invariant),
        "a5_invariant_conference_cubic_line_count": len(distinct_lines),
        "projective_line_match_count": len(matches),
        "matching_matrices": [record["matrix"] for record in matches],
        "projected_sheet_cubic_f11_singular_points": orbit_record(projected_singular, data),
        "projected_sheet_cubic_chordal_identification": chordal_certificate(
            data, data["projected_cubic"]
        ),
        "conference_triangle_cubic_f11_singular_points": orbit_record(conference_singular, data),
        "a5_invariant_cubic_dimension": 2,
        "outer_normalizer_pencil_action": outer_pencil_action(
            data, conference_line, data["projected_cubic"]
        ),
        "invariant_pencil_singular_census": pencil,
        "inputs": {path.name: file_record(path) for path in INPUT_PATHS},
        "scope": "exact finite-field reconnaissance; naturality and inverse are not proved here",
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write:
        OUTPUT.write_text(rendered)
        print(f"wrote {OUTPUT}")
    else:
        if not OUTPUT.exists() or OUTPUT.read_text() != rendered:
            raise SystemExit("stale certificate")
        print(f"CHECK OK ({result['verdict']})")


if __name__ == "__main__":
    main()
