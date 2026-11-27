# Changelog

## Unreleased

- `sign(p, n?)`, `isEven`, `isOdd` and `inversions` for a single
  permutation (array of images or cycle text).
- `constructible` decides polynomials of degree 8 and more when they are,
  coefficient for coefficient, the minimal polynomial of cos(2π/n) or
  2cos(2π/n) (Gauss–Wantzel); the factor carries `cosine: { n, scale }`.
- `constructiblePolygon(n)`, another name for `regularPolygon(n)`.

## 0.5.0

Presentations, abelian groups and modules (roadmap phase 5).

- `presentation("a, b | a^4, b^2, (ab)^2")`: groups from generators and
  relators or relations; words with powers, brackets, commutators and capital
  letters for inverses.
- Todd–Coxeter coset enumeration (HLT with coincidences): `P.enumerate(H)`,
  `P.index(H)`, `P.order()`, `P.toGroup()` (elements printed as shortest
  words). Stops with an error past `maxCosets`.
- `presentationOf(G)` for any finite Group; `abelianInvariants(G or P)`.
- `smithForm(A, domain)` over the integers (`ZZ`) or Q[x]
  (`polynomialDomain`), with U, V and V⁻¹.
- `AbelianGroup` from relations (`fromRelations`, `integerModule`), invariant
  factors (`abelianGroup`) or cyclic orders (`cyclicProduct`); elementary
  divisors, products, `toGroup()`, generators of the cyclic factors;
  `abelianGroups(n)`.
- `matrixModule(A)`: Qⁿ as a Q[x]-module, with invariant factors, elementary
  divisors, minimal and characteristic polynomials and the rational canonical
  form; `companionMatrix(f)`.
- Tests compare the Smith forms and module structures with @ratmath/linalg
  (now a submodule under vendor/, used by the tests only).
- Docs: pages on presentations, abelian groups and modules; the group
  explorer accepts presentations.

## 0.4.0

Fields and large groups (roadmap phase 6) and character tables (phase 7).

- `schreierSims(generators)` and `PermChain`: permutation groups held by a
  base and strong generating set (deterministic Schreier–Sims), with order,
  membership, words in the generators, random elements, orbits, stabilisers,
  subgroup comparisons, normal closures, the derived series, block systems
  and primitivity. `largeGroups` has Sₙ, Aₙ, the Mathieu groups and the
  Rubik's cube and pocket cube groups.
- `extension({ … })`: towers of field extensions of Q with exact elements,
  degrees and the tower law, minimal polynomials over Q or a level of the
  tower, norm and trace, primitive elements, roots and factors over the
  field, automorphisms, `galoisGroup()` and `fixedField()`.
- `solveByRadicals(poly)`: Galois group, solvability and, when solvable, the
  roots as nested radicals for irreducible factors up to degree 7, each
  formula checked against its root. `RadicalExpression` prints in steps, in
  one line or as TeX. `rootOfUnityByRadicals(n)`.
- `constructible(poly)`, `regularPolygon(n)` (Gauss–Wantzel, with cos(2π/n)
  in square roots) and `classicalProblems`.
- `conwayPolynomial(p, n)` and `conwayField(p, n)` with compatible
  embeddings.
- `characterTable(G)`: irreducible characters with exact values in a
  cyclotomic field, by the Burnside–Dixon–Schneider method (class matrices
  over GF(p), then Dixon's lift), checked by row and column orthogonality.
  `ClassFunction` has inner products, decomposition, products, symmetric and
  alternating squares, Adams operations, the Frobenius–Schur indicator,
  kernels, restriction, induction and permutation characters. `Cyclotomic`
  is the exact number type of the entries.
- Docs: two more topic pages (character tables, working with characters) and
  `characters` and `tensor` in the explorer.
- Docs: six new topic pages (large permutation groups, field extensions and
  towers, automorphisms and the Galois correspondence, solving by radicals,
  ruler and compass, Conway polynomials), and `large`, `tower`, `galois`,
  `radicals`, `unity`, `polygon`, `constructible` and `conway` in the
  explorer.

## 0.3.0

Rings, ideals and Gröbner bases (roadmap phases 3 and 4).

- `groebner(polys, { order, certificates, maxPairs, timeLimit })`: Buchberger's
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
