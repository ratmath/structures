// @ratmath/structures: abstract algebra for RatMath.
//
// Finite groups (named groups, subgroups, quotients, Sylow subgroups,
// homomorphisms, actions); polynomials in several variables with the
// rational functions built from them, over any exact field; Gröbner bases,
// ideals, quotient rings and polynomial systems; the integers mod n and
// quadratic integer rings.

export { ring, PolyRing, MPoly } from "./src/mpoly.js";
export { QQ } from "./src/field.js";
export { ORDERS } from "./src/monomial.js";
export { gcd, lcm, content, primitivePart, squareFree, squareFreePart } from "./src/gcd.js";
export { det, sylvester, resultant, discriminant, jacobian, hessian } from "./src/resultant.js";
export { factor, isIrreducible, integerPrimitive } from "./src/factor.js";
export { fractions, RationalFunction } from "./src/ratfun.js";
export { elementary, powerSum, completeHomogeneous, isSymmetric, toElementary } from "./src/symmetric.js";

// the fields polynomials are most often wanted over, from @ratmath/algebraic
export { primeField, galoisField, NumberField } from "@ratmath/algebraic";

// finite groups
export { Group, Coset } from "./src/group.js";
export {
  cyclic, dihedral, dicyclic, quaternion, kleinFour, symmetric, alternating, permutationGroup, fromPermGroup,
  matrixGroup, GL, SL, PSL, heisenberg, directProduct, semidirectProduct, metacyclic,
} from "./src/named.js";
export {
  Homomorphism, homomorphism, tryHomomorphism, identityMap, trivialMap, inclusion, projection,
  isomorphism, isIsomorphic, automorphisms, automorphismGroup, innerAutomorphism,
  firstIsomorphism, secondIsomorphism, thirdIsomorphism,
} from "./src/hom.js";
export { Action, action, regularAction, conjugationAction, cosetAction, burnside } from "./src/action.js";
export { smallGroups, numberOfGroups, identify, SMALL_ORDER_LIMIT } from "./src/small.js";
// algebraic's permutation groups (Galois groups are computed with them), re-exported
export { PermGroup, cycleString, fromCycles, cycleType } from "@ratmath/algebraic";

// Gröbner bases, ideals, quotient rings and polynomial systems
export { groebner, GroebnerBasis, sPolynomial, eliminationOrder } from "./src/groebner.js";
export { ideal, Ideal, chineseRemainder, QuotientRing, QuotientElement, solve } from "./src/ideal.js";

// the integers mod n and quadratic integer rings
export {
  crt, Zn, ZnRing, quadraticRing, QuadraticRing, QuadraticInteger, QuadraticIdeal, gaussianIntegers, eisensteinIntegers,
} from "./src/rings.js";
