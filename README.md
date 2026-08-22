# Chordal and Conference Cubics: Reconstruction and a Residual \(C_2\)-Torsor

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21895530-blue.svg)](https://doi.org/10.5281/zenodo.21895530)

[Read the paper (PDF).](chordal_conference_reconstruction.pdf)

Two geometrically inequivalent cubic shadows occur on the same
five-dimensional metric \(A_5\)-carrier: a chordal cubic and a conference
cubic. This paper asks whether they nevertheless retain equivalent source
information. It proves that they do after one chordal line is selected, and
computes the exact \(C_2\)-ambiguity left when that selection is forgotten.

The chordal cubic is singular along a rational normal quartic. Its twelve
\(\mathbf F_{11}\)-points form \(A_5/C_5\); the stabilizer quotient recovers
the original six-axis carrier \(A_5/D_{10}\). The normalized outer action
then gives an exact chosen-sign round trip between a selected chordal line and
the conference line. Without the selected line this is a residual double
cover, not an equivalence.

The recovered six-set also distinguishes the rank-five augmentation lattice
from the rank-six \(D\)-type weight lattice while identifying their common
binary heart. A uniform conference-saturation theorem explains the
split/inert mod-eight residue dichotomy. In order six the normalized conference
operator produces the unique nonsplit \(\mathbf F_4A_5\)-extension, and
global negation \(B\mapsto-B\) becomes Frobenius. Thus the
global-negation torsor of opposite conference switching classes—not the
distinct residual chordal-line torsor—is the Frobenius torsor.

This paper belongs to the *Clebsch: Rigidity from Sparse Shadows* series.
*Reconstructing the Clebsch Code from Its Deep-Hole Syndrome Locus*,
*Quadratic Trade Rigidity and Cubic Orientation in Conic Matching Quotients*,
and *Hitchin's Icosahedral Incidence Double Cover and Operator Realizations of
the Clebsch Cubic* supply concrete sources realizing the two cubic shadows.
*Reconstructing \(\operatorname{PG}(2,13)\), its conic, and polarity from the
minimum words of a binary conic code* is an independent minimum-word reconstruction branch whose
residual marking follows the same Frobenius-orbit principle.

## Building and checking

The manuscript is chordal_conference_reconstruction.tex. From this directory:

    make check

This replays the exact evidence program, builds the PDF, and fails on any
LaTeX warning. The checker uses only the Python standard library; the PDF
build uses a Nix-pinned TeX environment.

## Verification surface

The proof is human and structural. The program under
verification/evidence/ independently replays only the Paper-II normalization
leaf: the finite-field tensor placement and its forward/reverse identities.
The singular-locus, stabilizer, groupoid, lattice-saturation, extension, and
Frobenius arguments are printed in the paper and do not depend on the
program. See verification/README.md for the trust boundary. The adjacent
icosahedral-tower files are exploratory evidence for a separately gated
successor theorem and are not used by the manuscript check.

## License

CC BY 4.0. See LICENSE.
