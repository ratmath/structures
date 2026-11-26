# Changelog

## 0.2.0

Finite groups (roadmap phase 1).

- `Group`: any finite group held by its elements (a few thousand at most),
  from generators and a multiplication, a list of elements, or a Cayley table
  (axioms checked). Elements are given by value or by printed label.
- Named groups: `cyclic`, `dihedral` (order 2n), `dicyclic`, `quaternion`,
  `kleinFour`, `symmetric`, `alternating`, `permutationGroup`,
  `fromPermGroup`, `GL`, `SL`, `PSL`, `heisenberg`, `matrixGroup` over F_p;
  `directProduct`, `semidirectProduct` (the action is checked) and
  `metacyclic`. Each carries its natural action as `G.natural`.
- Inside a group: Cayley tables, element orders, exponent, conjugacy classes
  and the class equation, centre, centralisers, normalisers, commutators,
  derived, lower and upper central series, solvable, nilpotent, simple and
  perfect tests.
- Subgroups: every subgroup and the subgroup classes, normal subgroups (from
  conjugacy classes), intersections and joins, cosets, quotients,
  composition series, Sylow subgroups and their number.
- Homomorphisms from functions or images of generators (refused when the
  relations fail), kernels, images, preimages, composition; the three
  isomorphism theorems on examples; isomorphism tests; automorphism groups
  with their inner automorphisms.
- `smallGroups(n)` and `identify(G)`: the groups of order at most 20, one of
  each type, named.
- Actions: orbits, stabilisers, fixed points, kernels, the orbit-counting
  lemma (checked against the orbits), induced actions on subsets and
  colourings, regular, conjugation and coset actions, and `toPermGroup()` for
  Pólya counting with @ratmath/discrete.
- Re-exports algebraic's `PermGroup`, `cycleString`, `fromCycles`, `cycleType`.
- Docs: six new topic pages (groups, classes, subgroups, Sylow,
  homomorphisms, actions) and a group explorer on the overview, explore and
  topic pages; every group name on a page is checked by the build.

- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

First version: polynomials in several variables.

- Polynomial rings over Q or any RatMath field object, with lex, grlex,
  grevlex or custom monomial orders; reading and printing (Unicode, plain
  text, LaTeX).
- Arithmetic, derivatives, substitution, evaluation, composition,
  homogenisation.
- The division algorithm for one or several divisors.
- gcd, lcm, content, primitive part, square-free decomposition.
- Factoring over Q; irreducibility.
- Fraction-free determinants, Sylvester matrices, resultants, discriminants,
  Jacobians and Hessians.
- Rational functions as a field object (usable with @ratmath/linalg's field
  module for matrices with parameters).
- Symmetric polynomials and rewriting in elementary symmetric polynomials.
- Docs site with an overview, explorer, reference and seven topic pages.
