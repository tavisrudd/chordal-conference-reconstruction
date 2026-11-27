#!/usr/bin/env python3
"""Generate the paper-owned sparse-shadow export for Paper V."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import deque
from pathlib import Path


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"
NODE_PATH = EVIDENCE / "conference_node_completeness.json"
DET_PATH = EVIDENCE / "determinantal_presentation.json"
AXIS_PATH = EVIDENCE / "paper_ii_chordal_axis.json"
REQUIRED_EXPORT = "verification/evidence/sparse_shadow_export.json"


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def closure(generators: list[tuple[int, ...]], degree: int) -> set[tuple[int, ...]]:
    identity = tuple(range(degree))
    seen = {identity}
    pending = deque([identity])
    while pending:
        element = pending.popleft()
        for generator in generators:
            product = compose(generator, element)
            if product not in seen:
                seen.add(product)
                pending.append(product)
    return seen


def automorphisms(matrix: list[list[int]]) -> list[tuple[int, ...]]:
    degree = len(matrix)
    return [
        permutation
        for permutation in itertools.permutations(range(degree))
        if all(
            matrix[left][right] == matrix[permutation[left]][permutation[right]]
            for left in range(degree)
            for right in range(left + 1, degree)
        )
    ]


def rational(value: int) -> dict[str, int]:
    return {"denominator": 1, "numerator": value}


def compute() -> dict[str, object]:
    nodes = json.loads(NODE_PATH.read_text(encoding="utf-8"))
    determinant = json.loads(DET_PATH.read_text(encoding="utf-8"))
    axis = json.loads(AXIS_PATH.read_text(encoding="utf-8"))
    matrix = determinant["conference_matrix"]
    group = automorphisms(matrix)
    assert len(group) == 10
    action_generators = [
        (1, 0, 2, 3, 4, 5),
        (1, 2, 3, 4, 5, 0),
    ]
    assert len(closure(action_generators, 6)) == 720
    assert determinant["matrix_square_is_five_identity"]
    assert determinant["rank_at_most_one_locus_size"] == 6
    assert nodes["conference"]["f11_singular_point_count"] == 6
    assert nodes["chordal_control"]["f11_singular_point_count"] == 12
    assert nodes["conference"]["frame_points_are_the_singular_points"]
    outer = axis["outer_normalizer_pencil_action"]["outer_axis_permutation"]
    assert sorted(outer) == list(range(6))
    outer = tuple(outer)
    outer_square = compose(outer, outer)
    assert outer_square != tuple(range(6))
    assert compose(outer_square, outer_square) == tuple(range(6))
    assert axis["outer_normalizer_pencil_action"]["exchanged_chordal_parameters"] == [
        [0, 1],
        [1, 7],
    ]

    positive = []
    negative = []
    for left in range(6):
        for right in range(left + 1, 6):
            (positive if matrix[left][right] == 1 else negative).append([left, right])

    return {
        "schema": "sparse-shadow/v1",
        "profile": {
            "adapter": "paper_v_chordal_conference",
            "input": {
                "gate": {
                    "enabled": True,
                    "reason": "paper-owned sparse-shadow export frozen",
                    "required_export": REQUIRED_EXPORT,
                },
                "source": {
                    "paper": "V",
                    "theorem": "carrier recovery and exact marked return",
                    "artifact": "verification/evidence/conference_node_completeness.json",
                    "sha256": hashlib.sha256(NODE_PATH.read_bytes()).hexdigest(),
                },
                "base_field": {"kind": "rational"},
                "retained_residue": {
                    "action": "color_preserving_permutations",
                    "vertices": [
                        {"color": 0, "sign": 0, "weight": 1}
                        for _ in range(6)
                    ],
                    "relations": [
                        {"directed": False, "edges": positive, "name": "conference_positive"},
                        {"directed": False, "edges": negative, "name": "conference_negative"},
                    ],
                },
                "action": {
                    "kind": "vertex_permutations",
                    "degree": 6,
                    "generators": [list(value) for value in action_generators],
                },
                "selected_chordal_line": 0,
                "outer_involution": list(outer),
                "delta_matrix": [[rational(value) for value in row] for row in matrix],
                "verification_field": {
                    "characteristic": 11,
                    "degree": 1,
                    "modulus_coefficients_low_to_high": [],
                    "element_encoding": "least_nonnegative_residue",
                },
                "conference_singular_points": nodes["conference"]["f11_singular_points"],
                "chordal_singular_points": nodes["chordal_control"]["f11_singular_points"],
                "conference_cubic": nodes["conference"]["cubic"],
                "chordal_cubic": nodes["chordal_control"]["cubic"],
                "recovered_carrier": "singular quartic, twelve points, six axes, and conference switching class",
                "ambiguity": {"kind": "orientation_c2"},
                "odd_calibration": {
                    "name": "selected_chordal_line",
                    "support": [0],
                    "value": 1,
                },
                "minimality_collisions": [],
            },
        },
    }


def render(value: dict[str, object]) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", type=Path)
    mode.add_argument("--check", type=Path)
    args = parser.parse_args()
    content = render(compute())
    if args.write:
        args.write.write_text(content, encoding="utf-8")
    else:
        assert args.check.read_text(encoding="utf-8") == content
        print("Paper V sparse-shadow export: PASS")


if __name__ == "__main__":
    main()
