#!/usr/bin/env python3
"""Exact F_121 audit of the full stabilizer stratification on Paper V's quartic."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

P = 11
NON_SQUARE = 2
HERE = Path(__file__).resolve().parent
CORE_PATH = HERE / "paper_ii_chordal_axis.py"
OUTPUT_PATH = HERE / "icosahedral_stabilizer_tower.json"


def load_core():
    spec = importlib.util.spec_from_file_location("paper_v_core", CORE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {CORE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Elements of F_121 = F_11[u]/(u^2-2) are encoded as a + 11b.
def pair(x):
    return x % P, x // P


def elt(a, b=0):
    return (a % P) + P * (b % P)


def add(x, y):
    a, b = pair(x)
    c, d = pair(y)
    return elt(a + c, b + d)


def neg(x):
    a, b = pair(x)
    return elt(-a, -b)


def mul(x, y):
    a, b = pair(x)
    c, d = pair(y)
    return elt(a * c + NON_SQUARE * b * d, a * d + b * c)


def inv(x):
    a, b = pair(x)
    norm = (a * a - NON_SQUARE * b * b) % P
    if norm == 0:
        raise ZeroDivisionError
    scale = pow(norm, -1, P)
    return elt(a * scale, -b * scale)


def power(x, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, x)
        x = mul(x, x)
        exponent //= 2
    return result


def normalize(vector):
    pivot = next(value for value in vector if value)
    scale = inv(pivot)
    return tuple(mul(scale, value) for value in vector)


def base_matrix_vector(matrix, vector):
    return tuple(
        sum_field(mul(entry, value) for entry, value in zip(row, vector))
        for row in matrix
    )


def sum_field(values):
    result = 0
    for value in values:
        result = add(result, value)
    return result


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(right)))


def invert(permutation):
    result = [0] * len(permutation)
    for source, image in enumerate(permutation):
        result[image] = source
    return tuple(result)


def element_order(permutation):
    identity = tuple(range(len(permutation)))
    value = identity
    for order in range(1, 61):
        value = compose(permutation, value)
        if value == identity:
            return order
    raise AssertionError("group element order exceeds 60")


def build():
    core = load_core()
    data = core.h3_data()
    base_certificate = core.build_certificate()
    projectivity = base_certificate[
        "projected_sheet_cubic_chordal_identification"
    ]["projectivity"]
    group = tuple(data["parent_group"])
    actions = data["axis_actions"]

    parameters = [(1, t) for t in range(P * P)] + [(0, 1)]
    points = []
    for s, t in parameters:
        veronese = tuple(mul(power(s, 4 - exponent), power(t, exponent)) for exponent in range(5))
        points.append(normalize(base_matrix_vector(projectivity, veronese)))
    assert len(set(points)) == 122

    stabilizers = []
    for point in points:
        stabilizers.append(
            frozenset(
                g for g in group if normalize(base_matrix_vector(actions[g], point)) == point
            )
        )
    order_counts = {}
    for subgroup in stabilizers:
        order_counts[len(subgroup)] = order_counts.get(len(subgroup), 0) + 1
    assert order_counts == {5: 12, 3: 20, 2: 30, 1: 60}

    subgroup_records = {}
    normalizer_sizes = {}
    for subgroup in set(stabilizers):
        order = len(subgroup)
        subgroup_records.setdefault(order, []).append(subgroup)
        if order == 1:
            continue
        normalizer = [
            g
            for g in group
            if frozenset(
                compose(compose(g, h), invert(g)) for h in subgroup
            )
            == subgroup
        ]
        normalizer_sizes.setdefault(order, set()).add(len(normalizer))

    unique_subgroups = {order: len(groups) for order, groups in subgroup_records.items()}
    fibre_sizes = {
        order: sorted(stabilizers.count(subgroup) for subgroup in groups)
        for order, groups in subgroup_records.items()
    }
    assert unique_subgroups == {5: 6, 3: 10, 2: 15, 1: 1}
    assert fibre_sizes[5] == [2] * 6
    assert fibre_sizes[3] == [2] * 10
    assert fibre_sizes[2] == [2] * 15
    assert fibre_sizes[1] == [60]
    assert normalizer_sizes == {5: {10}, 3: {6}, 2: {4}}

    f11_indices = [
        index
        for index, point in enumerate(points)
        if all(pair(coordinate)[1] == 0 for coordinate in point)
    ]
    f11_stabilizer_counts = {}
    for index in f11_indices:
        order = len(stabilizers[index])
        f11_stabilizer_counts[order] = f11_stabilizer_counts.get(order, 0) + 1
    assert len(f11_indices) == 12
    assert f11_stabilizer_counts == {5: 12}

    free_indices = [index for index, subgroup in enumerate(stabilizers) if len(subgroup) == 1]
    first = points[free_indices[0]]
    free_orbit = {
        normalize(base_matrix_vector(actions[g], first))
        for g in group
    }
    assert len(free_orbit) == 60
    assert free_orbit == {points[index] for index in free_indices}

    class_orders = {}
    for g in group:
        order = element_order(g)
        class_orders[order] = class_orders.get(order, 0) + 1
    assert class_orders == {1: 1, 2: 15, 3: 20, 5: 24}

    return {
        "schema": "paper-v-icosahedral-stabilizer-tower-v1",
        "field_model": "F_11[u]/(u^2-2)",
        "f121_projective_line_points": len(points),
        "geometric_stabilizer_order_counts": {
            str(order): count for order, count in sorted(order_counts.items())
        },
        "unique_exact_stabilizer_subgroups": {
            str(order): count for order, count in sorted(unique_subgroups.items())
        },
        "fixed_pair_fibre_sizes": {
            str(order): sizes for order, sizes in sorted(fibre_sizes.items()) if order != 1
        },
        "normalizer_orders": {
            str(order): sorted(sizes) for order, sizes in sorted(normalizer_sizes.items())
        },
        "f11_rational_point_count": len(f11_indices),
        "f11_stabilizer_order_counts": {
            str(order): count for order, count in sorted(f11_stabilizer_counts.items())
        },
        "free_locus_is_one_regular_orbit_over_f121": True,
        "a5_element_order_counts": {
            str(order): count for order, count in sorted(class_orders.items())
        },
        "base_certificate_sha256": hashlib.sha256(
            (HERE / "paper_ii_chordal_axis.json").read_bytes()
        ).hexdigest(),
        "base_checker_sha256": hashlib.sha256(CORE_PATH.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    certificate = build()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check:
        assert OUTPUT_PATH.read_text() == rendered
        print("CHECK OK (12+20+30+60)")
    else:
        OUTPUT_PATH.write_text(rendered)
        print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
