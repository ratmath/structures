// @ratmath/structures: abstract algebra for RatMath.
//
// This first part is polynomials in several variables and the rational
// functions built from them, over any exact field.

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
