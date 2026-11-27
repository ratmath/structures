# Changelog

## 0.3.0

Rings, ideals and Gröbner bases (roadmap phases 3 and 4).

- `groebner(polys, { order, certificates, maxPairs })`: Buchberger's
  algorithm with the normal selection strategy and Gebauer and Möller's
  criteria, fraction-free over Q and by field arithmetic elsewhere. The
  result has the reduced basis, `reduce`, `divide`, `contains`,
  `express(f)` (a member of the ideal as a combination of the original
  generators, multiplied out to confirm), `verify()` (Buchberger's criterion
  re-checked), `dimension`, `standardMonomials()`, `quotientDimension` and
  `stats`. Also `sPolynomial` and `eliminationOrder`.
- `ideal(R, generators)`: membership, equality, `add`, `mul`, `pow`,
  `intersect`, `quotient`, `saturate`, `eliminate`, `radicalContains` (with
  a certificate fᵏ = Σ cᵢfᵢ on request), `exponentOf`, `radical`,
  `isRadical`, `dimension`, `degree`, `minimalPolynomial`, `isMaximal`,
  `isPrime`, and `chineseRemainder(I, J)`.
- `solve(R, polys, { complex })` and `I.solve()`: exact solutions of systems
  with finitely many solutions over Q, by a separating linear form and a
  rational univariate representation. Real solutions as algebraic numbers
  with float approximations, complex ones on request, multiplicity, and the
  components with their number fields.
- Quotient rings `I.quotientRing()`: elements in normal form with `inverse`,
  `isUnit`, `isZeroDivisor`, `isNilpotent`, `minimalPolynomial`, `matrix`;
  `isField`, `isIntegralDomain`, `isReduced`, `asField`.
- `Zn(n)`: arithmetic, units, zero divisors, nilpotents, idempotents, ideals,
  the unit group, primitive roots, `decompose()`; `crt(residues, moduli)`
  (moduli need not be coprime).
- `quadraticRing(d)`, `gaussianIntegers()`, `eisensteinIntegers()`: elements,
  units and the fundamental unit, `divmod` and `gcd`, `elementsOfNorm`,
  `isIrreducible` against `isPrimeElement`, `factorizations` (all of them),
  `factor` (when unique), ideals with `mul`, `isPrincipal`, `generator`,
  `factor`, `splitting(p)`, and `classGroup()`.
- Docs: six new topic pages (Gröbner bases, ideals, solving polynomial
  systems, quotient rings, the integers mod n, quadratic integer rings), and
  `groebner`, `member`, `eliminate`, `radical`, `solve`, `zmod`, `crt`,
  `quadratic` and `qfactor` in the explorer.

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
