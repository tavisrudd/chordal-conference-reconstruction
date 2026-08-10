#!/usr/bin/env python3
"""Exact certificate for the A5/D5 orbit code and q=11 syndrome scheme."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "a5_subgroup_decoder.json"

Vector = tuple[int, int, int]
Matrix = tuple[Vector, Vector, Vector]


@dataclass(frozen=True, order=True)
class QI:
    """a+b*tau in Z[tau], tau^2=tau+1."""

    a: int
    b: int = 0

    def __add__(self, other: QI | int) -> QI:
        other = other if isinstance(other, QI) else QI(other)
        return QI(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> QI:
        return QI(-self.a, -self.b)

    def __sub__(self, other: QI | int) -> QI:
        return self + (-other if isinstance(other, QI) else -QI(other))

    def __rsub__(self, other: QI | int) -> QI:
        return (-self) + other

    def __mul__(self, other: QI | int) -> QI:
        other = other if isinstance(other, QI) else QI(other)
        return QI(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    __rmul__ = __mul__

    def norm(self) -> int:
        return self.a * self.a + self.a * self.b - self.b * self.b

    def pair(self) -> list[int]:
        return [self.a, self.b]


def determinant(rows: list[list[QI]]) -> QI:
    total = QI(0)
    n = len(rows)
    for permutation in itertools.permutations(range(n)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(n)
            for j in range(i + 1, n)
        )
        term = QI(-1 if inversions % 2 else 1)
        for i, j in enumerate(permutation):
            term *= rows[i][j]
        total += term
    return total


def symbolic_points() -> list[tuple[QI, QI, QI]]:
    z, o, t = QI(0), QI(1), QI(0, 1)
    return [
        (z, o, o - t),
        (z, o, t - o),
        (o, o - t, z),
        (o, t - o, z),
        (o, z, -t),
        (o, z, t),
    ]


def det3_symbolic(columns: tuple[tuple[QI, QI, QI], ...]) -> QI:
    return determinant([[columns[j][i] for j in range(3)] for i in range(3)])


def veronese(point: tuple[QI, QI, QI]) -> list[QI]:
    x, y, z = point
    return [x * x, y * y, z * z, x * y, x * z, y * z]


def inv(value: int, q: int) -> int:
    assert value % q
    return pow(value, q - 2, q)


def normalize(vector: Vector, q: int) -> Vector:
    pivot = next(value % q for value in vector if value % q)
    scale = inv(pivot, q)
    return tuple(value * scale % q for value in vector)  # type: ignore[return-value]


def dot(left: Vector, right: Vector, q: int) -> int:
    return sum(a * b for a, b in zip(left, right)) % q


def mat_mul(left: Matrix, right: Matrix, q: int) -> Matrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3)) % q for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def mat_vec(matrix: Matrix, vector: Vector, q: int) -> Vector:
    return tuple(dot(row, vector, q) for row in matrix)  # type: ignore[return-value]


def mat_normalize(matrix: Matrix, q: int) -> Matrix:
    pivot = next(value % q for row in matrix for value in row if value % q)
    scale = inv(pivot, q)
    return tuple(tuple(value * scale % q for value in row) for row in matrix)  # type: ignore[return-value]


def mat_inverse(matrix: Matrix, q: int) -> Matrix:
    work = [list(row) + [int(i == j) for j in range(3)] for i, row in enumerate(matrix)]
    for column in range(3):
        pivot = next(row for row in range(column, 3) if work[row][column] % q)
        work[column], work[pivot] = work[pivot], work[column]
        scale = inv(work[column][column] % q, q)
        work[column] = [value * scale % q for value in work[column]]
        for row in range(3):
            if row != column:
                factor = work[row][column] % q
                work[row] = [
                    (left - factor * right) % q
                    for left, right in zip(work[row], work[column])
                ]
    return tuple(tuple(row[3:]) for row in work)  # type: ignore[return-value]


def columns_matrix(columns: list[Vector]) -> Matrix:
    return tuple(tuple(columns[j][i] for j in range(3)) for i in range(3))  # type: ignore[return-value]


def frame_matrix(points: list[Vector], q: int) -> Matrix:
    basis = columns_matrix(points[:3])
    coordinates = mat_vec(mat_inverse(basis, q), points[3], q)
    assert all(coordinates)
    return tuple(
        tuple(basis[i][j] * coordinates[j] % q for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def frame_map(source: list[Vector], target: list[Vector], q: int) -> Matrix:
    return mat_normalize(
        mat_mul(frame_matrix(target, q), mat_inverse(frame_matrix(source, q), q), q), q
    )


def h3_roots(q: int, tau: int) -> set[Vector]:
    roots = {(1, 0, 0), (0, 1, 0), (0, 0, 1)}
    for left_sign, right_sign in itertools.product((1, -1), repeat=2):
        root = (1, left_sign * tau % q, right_sign * (tau - 1) % q)
        roots.update(normalize(root[i:] + root[:i], q) for i in range(3))
    assert len(roots) == 15
    return roots


def six_points(q: int, tau: int) -> set[Vector]:
    raw = [
        (0, 1, 1 - tau),
        (0, 1, tau - 1),
        (1, 1 - tau, 0),
        (1, tau - 1, 0),
        (1, 0, -tau),
        (1, 0, tau),
    ]
    return {normalize(point, q) for point in raw}


def reflection(root: Vector, q: int) -> Matrix:
    denominator = dot(root, root, q)
    factor = 2 * inv(denominator, q) % q
    return tuple(
        tuple((int(i == j) - factor * root[i] * root[j]) % q for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def reflection_group(q: int, roots: set[Vector]) -> set[Matrix]:
    identity: Matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    generators = [mat_normalize(reflection(root, q), q) for root in sorted(roots)]
    group = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            child = mat_normalize(mat_mul(current, generator, q), q)
            if child not in group:
                group.add(child)
                queue.append(child)
                assert len(group) <= 60
    assert len(group) == 60
    return group


def projective_stabilizer_group(points: set[Vector], q: int) -> set[Matrix]:
    """Independent construction: test all 6! permutations of the arc."""
    ordered = sorted(points)
    group: set[Matrix] = set()
    for target in itertools.permutations(ordered):
        matrix = frame_map(ordered[:4], list(target[:4]), q)
        if all(normalize(mat_vec(matrix, ordered[i], q), q) == target[i] for i in range(6)):
            group.add(matrix)
    assert len(group) == 60
    return group


def matrix_order(matrix: Matrix, q: int) -> int:
    identity: Matrix = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    power = identity
    for order in range(1, 601):
        power = mat_mul(power, matrix, q)
        if power == identity:
            return order
    raise AssertionError("matrix order exceeds group order")


def krylov_orderings(points: set[Vector], q: int) -> int:
    count = 0
    for ordering in itertools.permutations(sorted(points)):
        matrix = frame_map(list(ordering[:4]), list(ordering[1:5]), q)
        if normalize(mat_vec(matrix, ordering[4], q), q) == ordering[5]:
            count += 1
    return count


def all_vectors(q: int) -> list[Vector]:
    return list(itertools.product(range(q), repeat=3))  # type: ignore[return-value]


def vector_orbits(group: set[Matrix], q: int) -> list[set[Vector]]:
    linear_group = {
        tuple(tuple(scale * value % q for value in row) for row in matrix)
        for matrix in group
        for scale in range(1, q)
    }
    assert len(linear_group) == 600
    unseen = set(all_vectors(q))
    orbits: list[set[Vector]] = []
    while unseen:
        seed = min(unseen)
        orbit = {mat_vec(matrix, seed, q) for matrix in linear_group}
        unseen -= orbit
        orbits.append(orbit)
    return orbits


def label_orbits(orbits: list[set[Vector]], roots: set[Vector], columns: set[Vector], q: int) -> list[tuple[str, set[Vector]]]:
    labelled: list[tuple[str, set[Vector]]] = []
    ordinary: list[set[Vector]] = []
    for orbit in orbits:
        representative = min(orbit)
        if representative == (0, 0, 0):
            labelled.append(("zero", orbit))
            continue
        point = normalize(representative, q)
        multiplicity = sum(dot(root, point, q) == 0 for root in roots)
        if point in columns:
            label = "column_D5"
        elif multiplicity == 0:
            label = "deep_hole_C5"
        elif multiplicity == 2:
            label = "double_V4"
        elif multiplicity == 3:
            label = "triple_S3"
        elif multiplicity == 1:
            ordinary.append(orbit)
            continue
        else:
            raise AssertionError((representative, multiplicity))
        labelled.append((label, orbit))
    for index, orbit in enumerate(sorted(ordinary, key=lambda item: min(item)), 1):
        labelled.append((f"single_secant_C2_{index}", orbit))
    order = {
        "zero": 0,
        "column_D5": 1,
        "triple_S3": 2,
        "deep_hole_C5": 3,
        "double_V4": 4,
        "single_secant_C2_1": 5,
        "single_secant_C2_2": 6,
        "single_secant_C2_3": 7,
    }
    return sorted(labelled, key=lambda item: order[item[0]])


def intersection_tensor(classes: list[set[Vector]], q: int) -> list[list[list[int]]]:
    class_of = {vector: index for index, orbit in enumerate(classes) for vector in orbit}
    tensor: list[list[list[int]]] = []
    for target in classes:
        matrices = []
        for left in classes:
            row = []
            for right_index in range(len(classes)):
                values = {
                    sum(class_of[tuple((w[d] - u[d]) % q for d in range(3))] == right_index for u in left)
                    for w in target
                }
                assert len(values) == 1
                row.append(values.pop())
            matrices.append(row)
        tensor.append(matrices)
    return tensor


def check_tensor(tensor: list[list[list[int]]], valencies: list[int]) -> None:
    rank = len(valencies)
    for i in range(rank):
        for j in range(rank):
            assert sum(tensor[k][i][j] * valencies[k] for k in range(rank)) == valencies[i] * valencies[j]
            for k in range(rank):
                assert tensor[k][i][j] == tensor[k][j][i]
                assert tensor[k][0][j] == int(k == j)
    # Associativity of A_i A_j = sum_k p_ij^k A_k.
    for i, j, ell, m in itertools.product(range(rank), repeat=4):
        left = sum(tensor[k][i][j] * tensor[m][k][ell] for k in range(rank))
        right = sum(tensor[k][j][ell] * tensor[m][i][k] for k in range(rank))
        assert left == right


def q9_add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] + right[0]) % 3, (left[1] + right[1]) % 3)


def q9_mul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return (
        (left[0] * right[0] + left[1] * right[1]) % 3,
        (left[0] * right[1] + left[1] * right[0] + left[1] * right[1]) % 3,
    )


def q9_pow(value: tuple[int, int], exponent: int) -> tuple[int, int]:
    answer = (1, 0)
    while exponent:
        if exponent & 1:
            answer = q9_mul(answer, value)
        value = q9_mul(value, value)
        exponent //= 2
    return answer


def q9_normalize(vector: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    pivot = next(value for value in vector if value != (0, 0))
    inverse = q9_pow(pivot, 7)
    return tuple(q9_mul(inverse, value) for value in vector)


def nonsplit_check() -> dict[str, object]:
    points = symbolic_points()
    reduced = {
        q9_normalize(tuple((coordinate.a % 3, coordinate.b % 3) for coordinate in point))
        for point in points
    }
    frobenius = {
        q9_normalize(tuple(q9_pow(coordinate, 3) for coordinate in point)) for point in reduced
    }
    minors = [det3_symbolic(tuple(choice)) for choice in itertools.combinations(points, 3)]
    return {
        "field": "F_9=F_3[tau]/(tau^2-tau-1)",
        "arc_minor_nonzero_count": sum((value.a % 3, value.b % 3) != (0, 0) for value in minors),
        "frobenius_preserves_six_set": frobenius == reduced,
    }


def certificate() -> dict[str, object]:
    points = symbolic_points()
    quadratic_determinant = determinant([veronese(point) for point in points])
    expected = QI(-64, 48)  # 16(3*tau-4)
    assert quadratic_determinant == expected
    minors = [det3_symbolic(tuple(choice)) for choice in itertools.combinations(points, 3)]
    assert len(minors) == 20 and all(value.norm() in {4, -4} for value in minors)

    q, tau = 11, 8
    roots = h3_roots(q, tau)
    columns = six_points(q, tau)
    reflections = reflection_group(q, roots)
    stabilizer = projective_stabilizer_group(columns, q)
    assert reflections == stabilizer
    assert all({normalize(mat_vec(matrix, point, q), q) for point in columns} == columns for matrix in reflections)

    raw_orbits = vector_orbits(reflections, q)
    labelled = label_orbits(raw_orbits, roots, columns, q)
    labels = [label for label, _ in labelled]
    classes = [orbit for _, orbit in labelled]
    assert [len(orbit) for orbit in classes] == [1, 60, 100, 120, 150, 300, 300, 300]
    tensor = intersection_tensor(classes, q)
    check_tensor(tensor, [len(orbit) for orbit in classes])

    linear_group = {
        tuple(tuple(scale * value % q for value in row) for row in matrix)
        for matrix in reflections
        for scale in range(1, q)
    }
    stabilizers = {}
    for label, orbit in labelled[1:]:
        representative = min(orbit)
        fixed = {matrix for matrix in linear_group if mat_vec(matrix, representative, q) == representative}
        stabilizers[label] = {
            "order": len(fixed),
            "element_orders": dict(sorted(Counter(matrix_order(matrix, q) for matrix in fixed).items())),
        }

    ordered_columns = sorted(columns)
    permutations = {
        tuple(ordered_columns.index(normalize(mat_vec(matrix, point, q), q)) for point in ordered_columns)
        for matrix in reflections
    }
    three_subsets = set(itertools.combinations(range(6), 3))
    subset_orbits = []
    while three_subsets:
        seed = min(three_subsets)
        orbit = {
            tuple(sorted(permutation[index] for index in seed)) for permutation in permutations
        }
        three_subsets -= orbit
        subset_orbits.append(sorted(orbit))
    assert sorted(map(len, subset_orbits)) == [10, 10]

    return {
        "schema": "a5_data-a5-subgroup-decoder-v1",
        "symbolic_orbit_mds": {
            "ring": "Z[tau]/(tau^2-tau-1)",
            "three_by_three_minor_norms": dict(sorted(Counter(value.norm() for value in minors).items())),
            "quadratic_evaluation_determinant": quadratic_determinant.pair(),
            "quadratic_evaluation_determinant_text": "16(3*tau-4)",
            "norm": quadratic_determinant.norm(),
            "bad_characteristics": [2, 5],
            "nonsplit_sample": nonsplit_check(),
        },
        "q11": {
            "projective_group_order": len(reflections),
            "independent_group_constructions_equal": True,
            "krylov_orderings_checked": 720,
            "krylov_orderings_passing": krylov_orderings(columns, q),
            "three_subset_A5_orbit_sizes": sorted(map(len, subset_orbits)),
            "syndrome_relation_labels": labels,
            "syndrome_relation_valencies": [len(orbit) for orbit in classes],
            "point_stabilizers": stabilizers,
            "intersection_matrices_by_target_relation": tensor,
            "intersection_algebra_commutative_associative": True,
        },
    }


def canonical_bytes(data: dict[str, object]) -> bytes:
    return (json.dumps(data, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = canonical_bytes(certificate())
    if args.write:
        OUTPUT.write_bytes(data)
    if args.check:
        assert OUTPUT.read_bytes() == data
    if not args.write and not args.check:
        print(data.decode(), end="")
    print(f"sha256={hashlib.sha256(data).hexdigest()} bytes={len(data)}")


if __name__ == "__main__":
    main()
