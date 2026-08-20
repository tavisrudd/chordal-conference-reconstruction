#!/usr/bin/env python3
"""Exact mod-11 completeness certificate for the conference node count.

Proposition (Nodes persist in characteristic eleven) proves structurally that the
six points p_a = [1 - 6 e_a] are ordinary nodes of the conference triangle cubic
c_B over the algebraic closure of F_11.  This program certifies the complementary
statement that there is nothing else: the singular subscheme of {c_B = 0} in
P(A_0) = P^4 over F_11 is finite of degree six, so the six reduced nodes exhaust
it after every base change.

The delegated step is a Groebner dimension-and-degree computation.  The reduction
that consumes it is human:

    projective degree 6, six distinct singular points, and a rank-four Hessian at
    each (hence local multiplicity one) force the singular scheme to be exactly
    those six reduced points over the algebraic closure.

The same computation is run on the chordal sheet cubic as a control; there the
singular scheme has Hilbert polynomial 4d + 1, the Hilbert polynomial of a
rational normal quartic curve, which carries no isolated point.

Everything is exact arithmetic over F_11 in the Python standard library.  The
conference and chordal cubics are read from the tracked certificate of
paper_ii_chordal_axis.py, and the conference cubic is additionally rebuilt here
from a pentagon-normalized order-six conference matrix, independently of that
program's projection route.

    python3 conference_node_completeness.py --check
    python3 conference_node_completeness.py --write
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
AXIS_CERTIFICATE = HERE / "paper_ii_chordal_axis.json"
OUTPUT = HERE / "conference_node_completeness.json"
SCHEMA = "paper-v-conference-node-completeness-v1"
P = 11
N = 5
BASIS = list(itertools.combinations_with_replacement(range(N), 3))
BASIS_INDEX = {monomial: index for index, monomial in enumerate(BASIS)}
INDEPENDENT_HILBERT_DEGREES = (3, 4, 5, 6)


# --------------------------------------------------------------------------
# the two cubics
# --------------------------------------------------------------------------
def pentagon_conference_matrix():
    """Symmetric order-six conference matrix, first row positive, pentagon negative."""
    matrix = [[0] * 6 for _ in range(6)]
    for index in range(1, 6):
        matrix[0][index] = matrix[index][0] = 1
    pentagon = {(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)}
    for left, right in itertools.combinations(range(1, 6), 2):
        sign = -1 if (left, right) in pentagon else 1
        matrix[left][right] = matrix[right][left] = sign
    return matrix


def matrix_square(matrix):
    size = len(matrix)
    return [
        [sum(matrix[i][k] * matrix[k][j] for k in range(size)) for j in range(size)]
        for i in range(size)
    ]


def triangle_cubic(matrix):
    """Restriction of the triangle cubic of `matrix` to the augmentation subspace.

    The augmentation is parametrized by y_a = x_a for a < 5 and x_5 = -(y_0+...+y_4);
    the returned vector lists coefficients against BASIS.
    """
    augmentation = [
        [(1 if axis == coordinate else 0) - (1 if axis == 5 else 0) for coordinate in range(N)]
        for axis in range(6)
    ]
    coefficients = [0] * len(BASIS)
    for i, j, k in itertools.combinations(range(6), 3):
        tau = matrix[i][j] * matrix[j][k] * matrix[k][i]
        if tau == 0:
            continue
        for a, b, c in itertools.product(range(N), repeat=3):
            index = BASIS_INDEX[tuple(sorted((a, b, c)))]
            coefficients[index] = (
                coefficients[index]
                + tau * augmentation[i][a] * augmentation[j][b] * augmentation[k][c]
            ) % P
    return coefficients


def projective_normalize(vector):
    for value in vector:
        if value % P:
            inverse = pow(value % P, P - 2, P)
            return [entry * inverse % P for entry in vector]
    return [entry % P for entry in vector]


def relabelling_to(target):
    """Axis permutation carrying the pentagon conference cubic onto `target`, if any."""
    base = pentagon_conference_matrix()
    for permutation in itertools.permutations(range(6)):
        relabelled = [[base[permutation[i]][permutation[j]] for j in range(6)] for i in range(6)]
        if projective_normalize(triangle_cubic(relabelled)) == target:
            return list(permutation)
    return None


# --------------------------------------------------------------------------
# polynomials over F_11 in five variables, degree-reverse-lexicographic
# --------------------------------------------------------------------------
def cubic_polynomial(coefficients):
    polynomial = {}
    for coefficient, monomial in zip(coefficients, BASIS):
        if coefficient % P:
            exponent = [0] * N
            for index in monomial:
                exponent[index] += 1
            polynomial[tuple(exponent)] = coefficient % P
    return polynomial


def total_degree(exponent):
    return sum(exponent)


def order_key(exponent):
    return (total_degree(exponent), tuple(-entry for entry in reversed(exponent)))


def lead(polynomial):
    return max(polynomial, key=order_key)


def exponent_lcm(left, right):
    return tuple(max(a, b) for a, b in zip(left, right))


def exponent_difference(left, right):
    return tuple(a - b for a, b in zip(left, right))


def divides(divisor, monomial):
    return all(a <= b for a, b in zip(divisor, monomial))


def shifted(polynomial, exponent, scalar):
    return {
        tuple(a + b for a, b in zip(monomial, exponent)): coefficient * scalar % P
        for monomial, coefficient in polynomial.items()
    }


def polynomial_sum(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        total = (result.get(monomial, 0) + coefficient) % P
        if total:
            result[monomial] = total
        elif monomial in result:
            del result[monomial]
    return result


def normal_form(polynomial, basis):
    polynomial = dict(polynomial)
    remainder = {}
    while polynomial:
        monomial = lead(polynomial)
        coefficient = polynomial[monomial]
        for element in basis:
            head = lead(element)
            if divides(head, monomial):
                factor = coefficient * pow(element[head], P - 2, P) % P
                polynomial = polynomial_sum(
                    polynomial,
                    shifted(element, exponent_difference(monomial, head), (-factor) % P),
                )
                break
        else:
            remainder[monomial] = coefficient
            del polynomial[monomial]
    return remainder


def s_polynomial(left, right):
    head_left, head_right = lead(left), lead(right)
    joint = exponent_lcm(head_left, head_right)
    scaled_left = shifted(
        left, exponent_difference(joint, head_left), pow(left[head_left], P - 2, P)
    )
    scaled_right = shifted(
        right, exponent_difference(joint, head_right), pow(right[head_right], P - 2, P)
    )
    return polynomial_sum(
        scaled_left, {monomial: (-value) % P for monomial, value in scaled_right.items()}
    )


def groebner_basis(generators):
    """Buchberger with the normal selection strategy and the coprime-lead criterion."""
    basis = [dict(generator) for generator in generators if generator]
    pairs = {(j, i) for i in range(len(basis)) for j in range(i)}
    while pairs:
        chosen = min(
            pairs,
            key=lambda pair: (
                total_degree(exponent_lcm(lead(basis[pair[0]]), lead(basis[pair[1]]))),
                exponent_lcm(lead(basis[pair[0]]), lead(basis[pair[1]])),
                pair,
            ),
        )
        pairs.discard(chosen)
        left, right = basis[chosen[0]], basis[chosen[1]]
        if all(min(a, b) == 0 for a, b in zip(lead(left), lead(right))):
            continue
        residue = normal_form(s_polynomial(left, right), basis)
        if residue:
            basis.append(residue)
            new = len(basis) - 1
            pairs.update((index, new) for index in range(new))
    return basis


def reduced_groebner_basis(basis):
    ordered = sorted(basis, key=lambda element: order_key(lead(element)))
    minimal = []
    for element in ordered:
        if any(divides(lead(kept), lead(element)) for kept in minimal):
            continue
        minimal.append(element)
    reduced = []
    for index, element in enumerate(minimal):
        others = minimal[:index] + minimal[index + 1:]
        residue = normal_form(element, others)
        scalar = pow(residue[lead(residue)], P - 2, P)
        reduced.append({m: c * scalar % P for m, c in residue.items()})
    return sorted(reduced, key=lambda element: order_key(lead(element)))


def verify_groebner(basis, generators):
    """Every generator and every S-pair reduces to zero against the basis."""
    for generator in generators:
        if normal_form(generator, basis):
            return False
    for left, right in itertools.combinations(basis, 2):
        if all(min(a, b) == 0 for a, b in zip(lead(left), lead(right))):
            continue
        if normal_form(s_polynomial(left, right), basis):
            return False
    return True


# --------------------------------------------------------------------------
# Hilbert series of the initial monomial ideal
# --------------------------------------------------------------------------
def minimal_monomials(monomials):
    ordered = sorted(set(monomials), key=total_degree)
    minimal = []
    for monomial in ordered:
        if not any(divides(kept, monomial) for kept in minimal):
            minimal.append(monomial)
    return tuple(minimal)


@functools.lru_cache(maxsize=None)
def hilbert_numerator(generators):
    """Numerator of the Hilbert series of S/(generators) over (1-t)^5."""
    if not generators:
        return (1,)
    if len(generators) == 1:
        degree = total_degree(generators[0])
        coefficients = [0] * (degree + 1)
        coefficients[0], coefficients[degree] = 1, -1
        return tuple(coefficients)
    pivot = generators[-1]
    rest = generators[:-1]
    head = hilbert_numerator(rest)
    quotient = minimal_monomials(
        exponent_difference(other, exponent_lcm_min(other, pivot)) for other in rest
    )
    tail = hilbert_numerator(quotient)
    degree = total_degree(pivot)
    result = [0] * max(len(head), len(tail) + degree)
    for index, value in enumerate(head):
        result[index] += value
    for index, value in enumerate(tail):
        result[index + degree] -= value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def exponent_lcm_min(left, right):
    return tuple(min(a, b) for a, b in zip(left, right))


def divide_by_one_minus_t(coefficients):
    quotient, carry = [], 0
    for value in coefficients[:-1]:
        carry += value
        quotient.append(carry)
    if coefficients[-1] + carry != 0:
        return None
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    return quotient


def hilbert_data(initial_generators):
    numerator = hilbert_numerator(minimal_monomials(initial_generators))
    reduced, pole = list(numerator), 0
    while True:
        step = divide_by_one_minus_t(reduced)
        if step is None:
            break
        reduced, pole = step, pole + 1
    codimension = pole
    projective_dimension = N - 1 - codimension
    degree = sum(reduced)
    return {
        "hilbert_numerator": list(numerator),
        "reduced_numerator": list(reduced),
        "one_minus_t_pole_order": codimension,
        "projective_dimension": projective_dimension,
        "projective_degree": degree,
        "hilbert_polynomial": hilbert_polynomial(reduced, N - codimension),
    }


def hilbert_polynomial(reduced, dimension):
    """Coefficients of HP(d) in the binomial-free monomial basis d^0, d^1, ...

    HP(d) = sum_i reduced[i] * C(d - i + dimension - 1, dimension - 1).
    Returned as the list of values HP(d) for d in a stated window plus the
    leading data, which is all the certificate needs.
    """
    if dimension == 0:
        return {"kind": "zero", "values": []}
    window = list(range(len(reduced) + dimension, len(reduced) + dimension + 6))
    values = []
    for d in window:
        total = 0
        for index, value in enumerate(reduced):
            total += value * binomial(d - index + dimension - 1, dimension - 1)
        values.append(total)
    return {"kind": "eventual", "degrees": window, "values": values}


def binomial(n, k):
    if k < 0 or n < 0 or n < k:
        return 0
    result = 1
    for step in range(k):
        result = result * (n - step) // (step + 1)
    return result


def hilbert_function_from_staircase(initial_generators, degree):
    generators = minimal_monomials(initial_generators)
    count = 0
    for monomial in monomials_of_degree(degree):
        if not any(divides(generator, monomial) for generator in generators):
            count += 1
    return count


def monomials_of_degree(degree):
    def walk(remaining, position):
        if position == N - 1:
            yield (remaining,)
            return
        for value in range(remaining + 1):
            for tail in walk(remaining - value, position + 1):
                yield (value,) + tail

    return walk(degree, 0)


# --------------------------------------------------------------------------
# independent Hilbert function: rank of the graded piece of the ideal
# --------------------------------------------------------------------------
def graded_hilbert_function(quadrics, degree):
    """dim (S/J)_degree computed by linear algebra alone, with no Groebner input."""
    monomials = list(monomials_of_degree(degree))
    position = {monomial: index for index, monomial in enumerate(monomials)}
    pivots = {}
    for quadric in quadrics:
        for shift in monomials_of_degree(degree - 2):
            row = [0] * len(monomials)
            for monomial, coefficient in quadric.items():
                row[position[tuple(a + b for a, b in zip(monomial, shift))]] = coefficient
            reduce_row(row, pivots)
    return len(monomials) - len(pivots)


def reduce_row(row, pivots):
    for column, pivot_row in pivots.items():
        if row[column]:
            factor = row[column]
            for index in range(len(row)):
                row[index] = (row[index] - factor * pivot_row[index]) % P
    for column, value in enumerate(row):
        if value:
            inverse = pow(value, P - 2, P)
            pivots[column] = [entry * inverse % P for entry in row]
            return


# --------------------------------------------------------------------------
# singular points and their Hessians
# --------------------------------------------------------------------------
def partial_derivative(polynomial, variable):
    result = {}
    for exponent, coefficient in polynomial.items():
        if exponent[variable]:
            lowered = list(exponent)
            lowered[variable] -= 1
            value = coefficient * exponent[variable] % P
            if value:
                key = tuple(lowered)
                result[key] = (result.get(key, 0) + value) % P
    return {key: value for key, value in result.items() if value}


def gradient(coefficients, point):
    values = [0] * N
    for coefficient, monomial in zip(coefficients, BASIS):
        if coefficient % P == 0:
            continue
        for position, index in enumerate(monomial):
            product = coefficient
            for other_position, other_index in enumerate(monomial):
                if other_position != position:
                    product = product * point[other_index] % P
            values[index] = (values[index] + product) % P
    return values


def hessian(coefficients, point):
    matrix = [[0] * N for _ in range(N)]
    for coefficient, monomial in zip(coefficients, BASIS):
        if coefficient % P == 0:
            continue
        for first, second in itertools.permutations(range(3), 2):
            row, column = monomial[first], monomial[second]
            remaining = monomial[3 - first - second]
            matrix[row][column] = (matrix[row][column] + coefficient * point[remaining]) % P
    return matrix


def matrix_rank(matrix):
    rows = [row[:] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((index for index in range(rank, len(rows)) if rows[index][column] % P), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], P - 2, P)
        rows[rank] = [entry * inverse % P for entry in rows[rank]]
        for index in range(len(rows)):
            if index != rank and rows[index][column] % P:
                factor = rows[index][column]
                rows[index] = [
                    (a - factor * b) % P for a, b in zip(rows[index], rows[rank])
                ]
        rank += 1
    return rank


def rational_singular_points(coefficients):
    points = []
    for pivot in range(N):
        for tail in itertools.product(range(P), repeat=N - 1 - pivot):
            point = (0,) * pivot + (1,) + tail
            if not any(gradient(coefficients, point)):
                points.append(point)
    return points


def frame_points():
    """The six points p_a = [1 - 6 e_a] of P(A_0) in augmentation coordinates."""
    points = []
    for axis in range(6):
        ambient = [(1 - 6 * (index == axis)) % P for index in range(6)]
        points.append(tuple(projective_normalize(ambient[:N])))
    return points


# --------------------------------------------------------------------------
# certificate
# --------------------------------------------------------------------------
def file_record(path):
    payload = path.read_bytes()
    return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def analyse(coefficients, independent_degrees):
    cubic = cubic_polynomial(coefficients)
    quadrics = [partial_derivative(cubic, variable) for variable in range(N)]
    basis = reduced_groebner_basis(groebner_basis(quadrics))
    initial = [lead(element) for element in basis]
    record = {
        "cubic": list(coefficients),
        "jacobian_generator_count": len(quadrics),
        "groebner_basis_size": len(basis),
        "groebner_basis_verified": verify_groebner(basis, quadrics),
        "initial_ideal_generators": [list(monomial) for monomial in minimal_monomials(initial)],
    }
    record.update(hilbert_data(initial))
    record["hilbert_function_staircase"] = {
        str(degree): hilbert_function_from_staircase(initial, degree)
        for degree in independent_degrees
    }
    record["hilbert_function_linear_algebra"] = {
        str(degree): graded_hilbert_function(quadrics, degree) for degree in independent_degrees
    }
    record["hilbert_functions_agree"] = (
        record["hilbert_function_staircase"] == record["hilbert_function_linear_algebra"]
    )
    points = rational_singular_points(coefficients)
    record["f11_singular_points"] = [list(point) for point in points]
    record["f11_singular_point_count"] = len(points)
    record["hessian_ranks"] = [matrix_rank(hessian(coefficients, point)) for point in points]
    record["ordinary_node_at_every_f11_singular_point"] = all(
        rank == N - 1 for rank in record["hessian_ranks"]
    )
    return record


def build_certificate():
    axis = json.loads(AXIS_CERTIFICATE.read_text())
    conference_cubic = axis["conference_triangle_cubic"]
    chordal_cubic = axis["projected_sheet_cubic"]

    matrix = pentagon_conference_matrix()
    square = matrix_square(matrix)
    permutation = relabelling_to(conference_cubic)

    conference = analyse(conference_cubic, INDEPENDENT_HILBERT_DEGREES)
    chordal = analyse(chordal_cubic, INDEPENDENT_HILBERT_DEGREES)

    frame = frame_points()
    conference["frame_points"] = [list(point) for point in frame]
    conference["frame_points_are_the_singular_points"] = set(frame) == {
        tuple(point) for point in conference["f11_singular_points"]
    }

    complete = (
        conference["groebner_basis_verified"]
        and conference["projective_dimension"] == 0
        and conference["projective_degree"] == 6
        and conference["f11_singular_point_count"] == 6
        and conference["ordinary_node_at_every_f11_singular_point"]
        and conference["frame_points_are_the_singular_points"]
        and conference["hilbert_functions_agree"]
        and chordal["projective_dimension"] == 1
        and chordal["projective_degree"] == 4
    )

    return {
        "schema": SCHEMA,
        "verdict": "EXACTLY_SIX_REDUCED_NODES" if complete else "INCOMPLETE",
        "field": P,
        "ambient": "P^4 = P(A_0) over F_11",
        "monomial_order": "degree-reverse-lexicographic",
        "independent_hilbert_degrees": list(INDEPENDENT_HILBERT_DEGREES),
        "conference_cubic_independent_construction": {
            "conference_matrix": matrix,
            "matrix_square_is_five_identity": square
            == [[5 if row == column else 0 for column in range(6)] for row in range(6)],
            "axis_permutation_onto_tracked_cubic": permutation,
            "reproduces_tracked_cubic": permutation is not None,
        },
        "conference": conference,
        "chordal_control": chordal,
        "inputs": {AXIS_CERTIFICATE.name: file_record(AXIS_CERTIFICATE)},
        "reduction": (
            "projective degree six, six distinct singular points, and a rank-four Hessian "
            "at each force the singular scheme of the conference cubic to be exactly those "
            "six reduced ordinary nodes over the algebraic closure of F_11"
        ),
        "scope": (
            "exact finite-field computation over F_11 only; it certifies no statement in "
            "characteristic zero and no statement about any other cubic in the pencil"
        ),
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
