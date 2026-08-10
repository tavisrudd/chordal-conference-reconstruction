# Verification

The manuscript's proof is human. The exact checker
evidence/paper_ii_chordal_axis.py independently replays the finite-field
linear algebra used as proof leaves:

- the five-isotypic projection of the Paper-II tensor;
- the outer-twisted projectivity and chordal Hankel equation;
- the scheme-support and twelve-point stabilizer census;
- recovery of the six original axes;
- the invariant cubic pencil and normalized outer action; and
- both forward/reverse tensor identities.

Run the release check from the repository root:

    make check

The accepted terminal line from the evidence program is CHECK OK. The
program uses exact arithmetic over \(\mathbf F_{11}\) and only the Python
standard library. The manuscript prints the conceptual arguments and marks
where this execution supplies exact matrix identities.

evidence/icosahedral_stabilizer_tower.py and its adjacent JSON file concern
an optional successor theorem about the \((2,3,5)\) inertia strata over
extension fields. They are retained as research evidence but are not invoked
by make check and certify no claim in this manuscript.
