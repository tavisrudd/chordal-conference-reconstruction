#!/usr/bin/env python3
"""Replay of the determinantal presentation of the conference cubic over F_11.

The node count of the conference triangle cubic is a theorem, proved without a
computer: the cubic is the determinant of the three-by-three matrix of linear
forms M(x) = sum_i x_i (l_i tensor m_i) built from the two eigenspaces of the
conference matrix B, and its rank-at-most-one locus is the six coordinate
tensors.  This program replays every step of that construction in the one
characteristic Paper V works in.  It is a cross-check of a proved statement, not
the evidence for it.

Checked here:

- B is a symmetric conference matrix and B^2 = 5I, with 5 a square in F_11;
- the eigenspaces W_+ and W_- are three-dimensional, mutually orthogonal, and
  are [6,3,4] maximum-distance-separable codes (every three-by-three minor
  nonzero, minimum weight four);
- sum_i l_i tensor m_i = 0, the relation that lets M descend to A_0;
- det M equals a nonzero multiple of the triangle cubic of B;
- M(1 - 6 e_a) = -6 l_a tensor m_a has rank one at each of the six frame points;
- the rank-at-most-one locus of M over P^4(F_11) is exactly those six points.

    python3 determinantal_presentation.py --check
    python3 determinantal_presentation.py --write
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import conference_node_completeness as base


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "determinantal_presentation.json"
SCHEMA = "paper-v-determinantal-presentation-v1"
P = base.P
N = base.N
BASIS = base.BASIS
BASIS_INDEX = base.BASIS_INDEX


def row_reduce(rows):
    rows = [row[:] for row in rows]
    pivots = {}
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column] % P), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], P - 2, P)
        rows[rank] = [entry * inverse % P for entry in rows[rank]]
        for index in range(len(rows)):
            if index != rank and rows[index][column] % P:
                factor = rows[index][column]
                rows[index] = [(a - factor * b) % P for a, b in zip(rows[index], rows[rank])]
        pivots[column] = rank
        rank += 1
    return rows, pivots, rank


def kernel_basis(matrix):
    """Canonical reduced-row-echelon kernel basis."""
    rows, pivots, _ = row_reduce(matrix)
    width = len(matrix[0])
    free = [column for column in range(width) if column not in pivots]
    basis = []
    for column in free:
        vector = [0] * width
        vector[column] = 1
        for pivot_column, pivot_row in pivots.items():
            vector[pivot_column] = (-rows[pivot_row][column]) % P
        basis.append(vector)
    return basis


def eigenspace(matrix, eigenvalue):
    shifted = [
        [(matrix[i][j] - (eigenvalue if i == j else 0)) % P for j in range(len(matrix))]
        for i in range(len(matrix))
    ]
    return kernel_basis(shifted)


def determinant3(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    ) % P


def rank_of(matrix):
    return row_reduce(matrix)[2]


def minimum_weight(basis):
    best = len(basis[0])
    for coefficients in itertools.product(range(P), repeat=len(basis)):
        if not any(coefficients):
            continue
        vector = [
            sum(coefficients[k] * basis[k][i] for k in range(len(basis))) % P
            for i in range(len(basis[0]))
        ]
        weight = sum(1 for entry in vector if entry)
        if weight < best:
            best = weight
    return best


def all_minors_nonzero(basis):
    return all(
        determinant3([[basis[r][c] for c in columns] for r in range(3)])
        for columns in itertools.combinations(range(6), 3)
    )


# --------------------------------------------------------------------------
# linear forms in the augmentation coordinates y_0..y_4, with x_5 = -(sum y)
# --------------------------------------------------------------------------
def ambient_to_augmentation(coefficients):
    """Rewrite a linear form sum_i a_i x_i on k^Omega in the coordinates y."""
    return [(coefficients[a] - coefficients[5]) % P for a in range(N)]


def multiply_linear(first, second):
    product = {}
    for i, a in enumerate(first):
        if a % P == 0:
            continue
        for j, b in enumerate(second):
            if b % P == 0:
                continue
            key = tuple(sorted((i, j)))
            product[key] = (product.get(key, 0) + a * b) % P
    return product


def multiply_quadratic_linear(quadratic, linear):
    cubic = [0] * len(BASIS)
    for (i, j), value in quadratic.items():
        if value % P == 0:
            continue
        for k, coefficient in enumerate(linear):
            if coefficient % P == 0:
                continue
            index = BASIS_INDEX[tuple(sorted((i, j, k)))]
            cubic[index] = (cubic[index] + value * coefficient) % P
    return cubic


def determinant_cubic(entries):
    """Expand det of a three-by-three matrix of linear forms in y."""
    cubic = [0] * len(BASIS)
    permutations = ((0, 1, 2), (1, 2, 0), (2, 0, 1), (0, 2, 1), (2, 1, 0), (1, 0, 2))
    signs = (1, 1, 1, -1, -1, -1)
    for permutation, sign in zip(permutations, signs):
        term = multiply_quadratic_linear(
            multiply_linear(entries[0][permutation[0]], entries[1][permutation[1]]),
            entries[2][permutation[2]],
        )
        cubic = [(a + sign * b) % P for a, b in zip(cubic, term)]
    return cubic


def build_certificate():
    matrix = base.pentagon_conference_matrix()
    square = base.matrix_square(matrix)
    reduced = [[value % P for value in row] for row in matrix]

    root = next(value for value in range(P) if value * value % P == 5 % P)
    plus = eigenspace(reduced, root)
    minus = eigenspace(reduced, (-root) % P)

    orthogonal = all(
        sum(u[i] * v[i] for i in range(6)) % P == 0 for u in plus for v in minus
    )

    # coordinate functionals in the bases above
    ell = [[plus[r][i] for r in range(3)] for i in range(6)]
    mee = [[minus[c][i] for c in range(3)] for i in range(6)]
    tensors = [[[ell[i][r] * mee[i][c] % P for c in range(3)] for r in range(3)] for i in range(6)]
    relation = [
        [sum(tensors[i][r][c] for i in range(6)) % P for c in range(3)] for r in range(3)
    ]

    # M(x) entrywise as linear forms in the augmentation coordinates
    entries = [
        [ambient_to_augmentation([tensors[i][r][c] for i in range(6)]) for c in range(3)]
        for r in range(3)
    ]
    det_cubic = determinant_cubic(entries)
    triangle = base.triangle_cubic(matrix)
    pivot = next((index for index, value in enumerate(det_cubic) if value % P), None)
    kappa = (
        det_cubic[pivot] * pow(triangle[pivot], P - 2, P) % P
        if pivot is not None and triangle[pivot] % P
        else None
    )
    proportional = kappa is not None and all(
        (kappa * a - b) % P == 0 for a, b in zip(triangle, det_cubic)
    )

    def matrix_at(y):
        ambient = list(y) + [(-sum(y)) % P]
        return [
            [sum(ambient[i] * tensors[i][r][c] for i in range(6)) % P for c in range(3)]
            for r in range(3)
        ]

    frame = base.frame_points()
    frame_records = []
    for axis, point in enumerate(frame):
        evaluated = matrix_at(point)
        scaled = [[(-6 * tensors[axis][r][c]) % P for c in range(3)] for r in range(3)]
        normalized = base.projective_normalize([entry for row in evaluated for entry in row])
        expected = base.projective_normalize([entry for row in scaled for entry in row])
        frame_records.append(
            {
                "axis": axis,
                "point": list(point),
                "rank": rank_of(evaluated),
                "equals_minus_six_times_coordinate_tensor": normalized == expected,
            }
        )

    low_rank = []
    for pivot_index in range(N):
        for tail in itertools.product(range(P), repeat=N - 1 - pivot_index):
            point = (0,) * pivot_index + (1,) + tail
            if rank_of(matrix_at(point)) <= 1:
                low_rank.append(list(point))

    mds = {
        "plus_all_three_by_three_minors_nonzero": all_minors_nonzero(plus),
        "minus_all_three_by_three_minors_nonzero": all_minors_nonzero(minus),
        "plus_minimum_weight": minimum_weight(plus),
        "minus_minimum_weight": minimum_weight(minus),
    }

    complete = (
        square == [[5 if r == c else 0 for c in range(6)] for r in range(6)]
        and len(plus) == 3
        and len(minus) == 3
        and orthogonal
        and relation == [[0] * 3 for _ in range(3)]
        and proportional
        and all(record["rank"] == 1 for record in frame_records)
        and all(record["equals_minus_six_times_coordinate_tensor"] for record in frame_records)
        and mds["plus_all_three_by_three_minors_nonzero"]
        and mds["minus_all_three_by_three_minors_nonzero"]
        and mds["plus_minimum_weight"] == 4
        and mds["minus_minimum_weight"] == 4
        and sorted(low_rank) == sorted(list(point) for point in frame)
    )

    return {
        "schema": SCHEMA,
        "verdict": "DETERMINANTAL_PRESENTATION_CONFIRMED" if complete else "INCOMPLETE",
        "field": P,
        "conference_matrix": matrix,
        "matrix_square_is_five_identity": square
        == [[5 if r == c else 0 for c in range(6)] for r in range(6)],
        "square_root_of_five": root,
        "eigenspace_plus_basis": plus,
        "eigenspace_minus_basis": minus,
        "eigenspaces_orthogonal": orthogonal,
        "maximum_distance_separable": mds,
        "coordinate_tensors": tensors,
        "coordinate_tensor_sum": relation,
        "determinant_cubic": det_cubic,
        "triangle_cubic": [value % P for value in triangle],
        "determinant_over_triangle_ratio": kappa,
        "determinant_equals_triangle_multiple": proportional,
        "frame_points": frame_records,
        "rank_at_most_one_locus": sorted(low_rank),
        "rank_at_most_one_locus_size": len(low_rank),
        "role": (
            "cross-check of the determinantal proof of the node count; the theorem is "
            "proved without a computer and holds in every characteristic outside "
            "{2,3,5} in which five is a square"
        ),
        "scope": "exact arithmetic over F_11 only; proves no statement on its own",
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
