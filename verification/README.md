# Verification

The manuscript's proof is human. The exact checker
evidence/paper_ii_chordal_axis.py independently replays the finite-field
linear algebra used in the one source-dependent normalization leaf:

- the five-isotypic projection of the Paper-II tensor;
- the outer-twisted projectivity and chordal Hankel equation;
- the invariant cubic pencil and normalized outer action; and
- both forward/reverse tensor identities.

Run the release check from the repository root:

    make check

The accepted terminal line from the evidence program is CHECK OK. The
program uses exact arithmetic over \(\mathbf F_{11}\) and only the Python
standard library. It is replay evidence, not a premise of the classification.
In particular, the manuscript proves the Hankel Jacobian saturation, the
special orbit and stabilizer quotient, the exact marking fibres, the uniform
root--weight saturation theorem, the nonsplit extension, and the Frobenius
identification without invoking this program.

The second checker evidence/conference_node_completeness.py certifies the
completeness half of the node count in characteristic eleven. The manuscript
proves structurally that the six frame points are ordinary nodes of the
conference cubic; the checker computes the Groebner dimension and degree of
the conference Jacobian ideal over \(\mathbf F_{11}\) and finds a
zero-dimensional singular scheme of degree six, with a rank-four Hessian at
each of the six \(\mathbf F_{11}\)-rational singular points. Degree six
against six multiplicity-one points is what forces exactness after base
change to \(\overline{\mathbf F}_{11}\); that reduction is human, and only
the dimension-and-degree step is delegated. The same computation on the
chordal sheet cubic returns Hilbert polynomial \(4d+1\), leaving no room for
an isolated point beside the rational normal quartic. This checker also uses
exact arithmetic over \(\mathbf F_{11}\) and only the Python standard
library, and it rebuilds the conference cubic from a pentagon-normalized
conference matrix rather than reusing the other program's projection route.
Its Hilbert function is computed twice, once from the initial-ideal staircase
and once by linear algebra on the graded pieces of the Jacobian ideal, with
no Groebner input. evidence/conference_node_completeness.sing is an
independent computer-algebra replay of the same two ideals; it needs Singular
and is therefore not invoked by make check.

The third checker evidence/determinantal_presentation.py replays the
determinantal presentation behind the node count. The conference cubic is the
determinant of the three-by-three matrix of linear forms
\(M(x)=\sum_i x_i\,\ell_i\otimes m_i\) assembled from the coordinate functionals
on the two eigenspaces of the conference matrix, and its rank-at-most-one locus
is the six coordinate tensors, which are the six nodes. That statement is a
theorem proved without a computer and valid in every characteristic outside
\(\{2,3,5\}\) in which five is a square; this checker replays it over
\(\mathbf F_{11}\), confirming the eigenspace splitting, the
maximum-distance-separable property of both eigenspaces, the relation
\(\sum_i\ell_i\otimes m_i=0\), the identity \(\det M=c_B\) for the basis it
fixes, and the rank-at-most-one locus. It is a cross-check of a proved
statement, not evidence for one, and the Groebner certificate above holds the
same status now that the count is proved structurally.

evidence/icosahedral_stabilizer_tower.py and its adjacent JSON file concern
an optional successor theorem about the \((2,3,5)\) inertia strata over
extension fields. They are retained as research evidence but are not invoked
by make check and certify no claim in this manuscript.
