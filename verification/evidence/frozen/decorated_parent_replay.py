#!/usr/bin/env python3
"""Independent replay of decorated-parent q=11 extension, iteration, and fibre claims."""

from __future__ import annotations

import itertools
import json
from collections import Counter, deque
from pathlib import Path

Q = 11
I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
J = ((1, 0, 0), (0, 0, 10), (0, 10, 0))
CERT = Path(__file__).with_name("decorated_parent.json")


def inv(value):
    return pow(value % Q, Q - 2, Q)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right)) % Q


def normalize(vector):
    pivot = next(value for value in vector if value % Q)
    scale = inv(pivot)
    return tuple(scale * value % Q for value in vector)


def mm(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) % Q for j in range(3))
        for i in range(3)
    )


def mv(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def normm(matrix):
    pivot = next(value for row in matrix for value in row if value % Q)
    scale = inv(pivot)
    return tuple(tuple(scale * value % Q for value in row) for row in matrix)


def closure(generators):
    answer = {I}
    queue = deque([I])
    while queue:
        left = queue.popleft()
        for right in generators:
            child = normm(mm(left, right))
            if child not in answer:
                answer.add(child)
                queue.append(child)
    return answer


def roots(tau):
    answer = {(1, 0, 0), (0, 1, 0), (0, 0, 1)}
    for sign1, sign2 in itertools.product((1, -1), repeat=2):
        vector = (1, sign1 * tau % Q, sign2 * (tau - 1) % Q)
        answer |= {normalize(vector[offset:] + vector[:offset]) for offset in range(3)}
    assert len(answer) == 15
    return answer


def reflection(vector):
    factor = 2 * inv(dot(vector, vector)) % Q
    return tuple(
        tuple((int(i == j) - factor * vector[i] * vector[j]) % Q for j in range(3))
        for i in range(3)
    )


def a5(tau):
    return closure([normm(reflection(vector)) for vector in sorted(roots(tau))])


def six_points(tau):
    return frozenset(
        normalize(vector)
        for vector in (
            (0, 1, 1 - tau),
            (0, 1, tau - 1),
            (1, 1 - tau, 0),
            (1, tau - 1, 0),
            (1, 0, -tau),
            (1, 0, tau),
        )
    )


def projective_points():
    return sorted(
        {
            normalize(vector)
            for vector in itertools.product(range(Q), repeat=3)
            if vector != (0, 0, 0)
        }
    )


def determinant(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    ) % Q


def cross(left, right):
    return normalize(
        (
            left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0],
        )
    )


def secants(points, plane):
    lines = {cross(left, right) for left, right in itertools.combinations(points, 2)}
    covered = {point for point in plane if any(dot(line, point) == 0 for line in lines)}
    return lines, covered


def rref(rows):
    matrix = [[value % Q for value in row] for row in rows]
    pivots = []
    row = 0
    if not matrix:
        return matrix, pivots
    for column in range(len(matrix[0])):
        choice = next((index for index in range(row, len(matrix)) if matrix[index][column]), None)
        if choice is None:
            continue
        matrix[row], matrix[choice] = matrix[choice], matrix[row]
        scale = inv(matrix[row][column])
        matrix[row] = [scale * value % Q for value in matrix[row]]
        for other in range(len(matrix)):
            if other != row and matrix[other][column]:
                multiple = matrix[other][column]
                matrix[other] = [
                    (value - multiple * pivot) % Q
                    for value, pivot in zip(matrix[other], matrix[row])
                ]
        pivots.append(column)
        row += 1
        if row == len(matrix):
            break
    return matrix[:row], pivots


def nullspace(rows, width):
    reduced, pivots = rref(rows)
    answer = []
    for free in range(width):
        if free in pivots:
            continue
        vector = [0] * width
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free] % Q
        answer.append(tuple(vector))
    return answer


def conic_row(point):
    x, y, z = point
    return (x * x % Q, y * y % Q, z * z % Q, x * y % Q, x * z % Q, y * z % Q)


def obstruction_matching(parent, child_conic):
    pairs = set()
    for omitted in parent:
        rows = [conic_row(point) for point in sorted(parent - {omitted})]
        kernel = nullspace(rows, 6)
        assert len(kernel) == 1
        coefficients = normalize(kernel[0])
        pair = frozenset(point for point in child_conic if dot(coefficients, conic_row(point)) == 0)
        assert len(pair) == 2
        pairs.add(pair)
    assert len(pairs) == 6 and set().union(*map(set, pairs)) == set(child_conic)
    return frozenset(pairs)


def matching_key(matching):
    return tuple(sorted(tuple(sorted(pair)) for pair in matching))


def encode_matching(matching):
    return [[list(point) for point in pair] for pair in matching_key(matching)]


def act_on_matching(matrix, matching):
    return frozenset(frozenset(normalize(mv(matrix, point)) for point in pair) for pair in matching)


def minimum_distance(generator):
    minimum = len(generator[0]) + 1
    count = 0
    for coefficients in itertools.product(range(Q), repeat=len(generator)):
        if not any(coefficients):
            continue
        word = [sum(coefficients[i] * generator[i][j] for i in range(len(generator))) % Q for j in range(len(generator[0]))]
        weight = sum(value != 0 for value in word)
        if weight < minimum:
            minimum, count = weight, 1
        elif weight == minimum:
            count += 1
    return minimum, count


def image(matrix, points):
    return frozenset(normalize(mv(matrix, point)) for point in points)


def conjugate(matrix, group):
    transpose = tuple(zip(*matrix))
    return frozenset(normm(mm(mm(matrix, element), transpose)) for element in group)


def commutator(left, right):
    left_inverse = tuple(zip(*left))
    right_inverse = tuple(zip(*right))
    return normm(mm(mm(mm(left, right), left_inverse), right_inverse))


def main():
    cert = json.loads(CERT.read_text())
    assert cert["schema"] == "othello.decorated.clebsch_deep_hole_extension.v2"
    plane = projective_points()
    plus = six_points(8)
    minus = six_points(4)
    conic = frozenset(point for point in plane if dot(point, point) == 0)
    assert [list(point) for point in sorted(plus)] == cert["tau8_parent"]
    assert [list(point) for point in sorted(minus)] == cert["tau4_parent"]
    assert [list(point) for point in sorted(conic)] == cert["deep_hole_conic"]
    assert frozenset(plane) - secants(plus, plane)[1] == conic

    replay_extensions = []
    for point in sorted(conic):
        seven = sorted(plus | {point})
        assert all(determinant(triple) for triple in itertools.combinations(seven, 3))
        check = [[point[column] for point in seven] for column in range(3)]
        kernel = nullspace(check, 7)
        assert len(kernel) == 4
        kernel_result = minimum_distance(kernel)
        dual_result = minimum_distance(check)
        bad = []
        for omitted in range(7):
            selected = [candidate for index, candidate in enumerate(seven) if index != omitted]
            if len(rref([conic_row(candidate) for candidate in selected])[1]) < 6:
                bad.append(list(seven[omitted]))
        assert len(bad) == 1 and bad[0] != list(point)
        replay_extensions.append((list(point), kernel_result, dual_result, bad))

    for replay, recorded in zip(replay_extensions, cert["extensions"]):
        point, kernel_result, dual_result, bad = replay
        assert point == recorded["point"]
        assert kernel_result == (recorded["kernel_parameters"][2], recorded["kernel_minimum_word_count"])
        assert dual_result == (recorded["dual_parameters"][2], recorded["dual_minimum_word_count"])
        assert bad == [item["omitted_point"] for item in recorded["bad_six_subsets"]]

    child_lines, child_covered = secants(conic, plane)
    assert (len(child_lines), len(child_covered), len(set(plane) - child_covered)) == (66, 133, 0)
    assert cert["second_transform"]["output_point_count"] == 0

    group = a5(8)
    assert len(group) == 60
    full = closure(list(group) + [J])
    assert len(full) == 1320 and all(image(matrix, conic) == conic for matrix in full)
    normalizer = [matrix for matrix in full if conjugate(matrix, group) == frozenset(group)]
    conjugates = {conjugate(matrix, group) for matrix in full}
    arcs = {image(matrix, plus) for matrix in full}
    assert (len(normalizer), len(conjugates), len(arcs)) == (60, 22, 22)
    assert all(frozenset(plane) - secants(arc, plane)[1] == conic for arc in arcs)

    matchings = {arc: obstruction_matching(arc, conic) for arc in arcs}
    assert len(set(matchings.values())) == 22
    plus_matching = matchings[plus]
    minus_matching = matchings[minus]
    matching_stabilizer = {
        matrix
        for matrix in full
        if act_on_matching(matrix, plus_matching) == plus_matching
    }
    assert matching_stabilizer == group
    assert act_on_matching(J, plus_matching) == minus_matching
    assert all(
        matchings[image(matrix, arc)] == act_on_matching(matrix, matching)
        for matrix in full
        for arc, matching in matchings.items()
    )
    assert cert["decorated_transform"]["distinct_obstruction_matching_count"] == 22
    assert cert["decorated_transform"]["injective_on_conjugate_a5_parent_locus"] is True
    assert cert["decorated_transform"]["pgl_equivariant"] is True

    # Reconstruct the index-two subgroup independently as the commutator subgroup of PGL_2(11).
    psl = closure([commutator(matrix, J) for matrix in full])
    assert len(psl) == 660 and psl < full and J not in psl
    remaining = set(arcs)
    sheets = []
    while remaining:
        representative = min(remaining, key=lambda arc: tuple(sorted(arc)))
        sheet = frozenset(image(matrix, representative) for matrix in psl)
        remaining -= sheet
        sheets.append(sheet)
    assert sorted(map(len, sheets)) == [11, 11]
    plus_sheet = next(sheet for sheet in sheets if plus in sheet)
    minus_sheet = next(sheet for sheet in sheets if minus in sheet)
    assert plus_sheet != minus_sheet
    assert frozenset(image(J, arc) for arc in plus_sheet) == minus_sheet

    plus_factorization = sorted((matchings[arc] for arc in plus_sheet), key=matching_key)
    minus_factorization = sorted((matchings[arc] for arc in minus_sheet), key=matching_key)
    all_edges = frozenset(frozenset(pair) for pair in itertools.combinations(conic, 2))
    for factorization in (plus_factorization, minus_factorization):
        multiplicities = Counter(edge for matching in factorization for edge in matching)
        assert frozenset(multiplicities) == all_edges
        assert Counter(multiplicities.values()) == Counter({1: 66})
    assert {act_on_matching(J, matching) for matching in plus_factorization} == set(minus_factorization)

    intersections = [[len(left & right) for right in minus_factorization] for left in plus_factorization]
    assert Counter(value for row in intersections for value in row) == Counter({1: 66, 0: 55})
    incidence = [[int(value == 1) for value in row] for row in intersections]
    columns = list(zip(*incidence))
    assert {sum(row) for row in incidence} == {6} and {sum(column) for column in columns} == {6}
    assert {
        sum(a * b for a, b in zip(incidence[i], incidence[j]))
        for i, j in itertools.combinations(range(11), 2)
    } == {3}
    assert {
        sum(a * b for a, b in zip(columns[i], columns[j]))
        for i, j in itertools.combinations(range(11), 2)
    } == {3}
    complement = [[1 - value for value in row] for row in incidence]
    complement_columns = list(zip(*complement))
    assert {sum(row) for row in complement} == {5}
    assert {sum(column) for column in complement_columns} == {5}
    assert {
        sum(a * b for a, b in zip(complement[i], complement[j]))
        for i, j in itertools.combinations(range(11), 2)
    } == {2}
    assert {
        sum(a * b for a, b in zip(complement_columns[i], complement_columns[j]))
        for i, j in itertools.combinations(range(11), 2)
    } == {2}

    recorded_factorization = cert["one_factorization_biplane"]
    assert recorded_factorization["psl_subgroup_order"] == 660
    assert recorded_factorization["pgl_over_psl_index"] == 2
    assert recorded_factorization["parent_sheet_sizes"] == [11, 11]
    assert recorded_factorization["matching_sheet_sizes"] == [11, 11]
    assert recorded_factorization["child_edge_count"] == 66
    assert recorded_factorization["edge_multiplicity_spectra"] == [{"1": 66}, {"1": 66}]
    assert recorded_factorization["tau8_sheet_matchings"] == [
        encode_matching(matching) for matching in plus_factorization
    ]
    assert recorded_factorization["tau4_sheet_matchings"] == [
        encode_matching(matching) for matching in minus_factorization
    ]
    assert recorded_factorization["golden_J_exchanges_sheets"] is True
    assert recorded_factorization["cross_matching_intersection_histogram"] == {"0": 55, "1": 66}
    assert recorded_factorization["share_edge_design_parameters"] == [11, 6, 3]
    assert recorded_factorization["disjointness_biplane_parameters"] == [11, 5, 2]
    assert recorded_factorization["share_edge_incidence_matrix"] == incidence

    unseen = set(arcs)
    orbit_records = []
    while unseen:
        representative = min(unseen, key=lambda arc: tuple(sorted(arc)))
        orbit = {image(matrix, representative) for matrix in group}
        unseen -= orbit
        orbit_records.append((len(orbit), plus in orbit, minus in orbit))
    orbit_records.sort()
    assert [item[0] for item in orbit_records] == [1, 5, 6, 10]
    assert next(size for size, has_plus, _ in orbit_records if has_plus) == 1
    assert next(size for size, _, has_minus in orbit_records if has_minus) == 5
    assert cert["marked_fibre"]["same_fibre_fixed_a5_orbit_sizes"] == [1, 5, 6, 10]
    assert cert["marked_fibre"]["golden_pair_is_complete"] is False
    print("independent replay: 12 MDS extensions; 12 A1 roots; D^2 empty; two K12 one-factorizations; biplane 2-(11,5,2)")


if __name__ == "__main__":
    main()
