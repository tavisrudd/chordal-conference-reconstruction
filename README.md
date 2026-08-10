# The Golden Companion Correspondence

[Read the paper (PDF).](golden_companion_reconstruction.pdf)

Paper II of the *Clebsch: Rigidity from Sparse Shadows* series produces a
signed cubic from a matching quotient. This paper identifies its canonical
five-dimensional residue and proves that it is a chordal companion—not the
conference cubic reconstructed in Papers I and III—in the same
two-dimensional \(A_5\)-invariant pencil.

The chordal cubic is singular along a rational normal quartic. Its twelve
\(\mathbf F_{11}\)-points form \(A_5/C_5\); pairing equal stabilizers recovers
the original six-axis carrier \(A_5/D_5\). The normalized outer action then
gives an exact oriented round trip between a selected chordal line and the
conference line. The theorem is marked and image-restricted: it retains the
labels and normalizations needed for literal source returns without claiming
that the two cubics are equal.

This is Paper V of the five-paper series. Papers I--III supply the marked
inputs, while Paper IV is an independent minimum-word reconstruction branch.

## Building and checking

The manuscript is golden_companion_reconstruction.tex. From this directory:

    make check

This replays the exact evidence program, builds the PDF, and fails on any
LaTeX warning. The checker uses only the Python standard library; the PDF
build uses a Nix-pinned TeX environment.

## Verification surface

The proof is human and structural. The program under
verification/evidence/ checks the exact finite-field tensor placement,
singular quartic, axis recovery, invariant pencil, normalized outer action,
and forward/reverse identities. See verification/README.md for the trust
boundary. The adjacent icosahedral-tower files are exploratory evidence for a
separately gated successor theorem and are not used by the manuscript check.

## License

CC BY 4.0. See LICENSE.
