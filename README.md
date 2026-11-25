# @ratmath/structures

Abstract algebra for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/). This first part is **polynomials
in several variables** with exact coefficients, and what is built directly on
them: division, greatest common divisors, factoring, resultants, rational
functions and symmetric polynomials.

Docs, with the mathematics explained and every example runnable:
<https://ratmath.com/structures/>

```js
import { ring, gcd, factor, resultant, discriminant, fractions, toElementary, powerSum } from "@ratmath/structures";

const R = ring("x, y, z");                       // Q[x, y, z]
const f = R.poly("x^3 + y^3 + z^3 - 3xyz");
factor(f).factors.map((t) => String(t.factor));  // ["x + y + z", "x² − xy + y² − xz − yz + z²"]
String(gcd(R.poly("x^3 - y^3"), R.poly("x^2 - y^2")));             // "x − y"
String(resultant(R.poly("x^2 + y^2 - 1"), R.poly("y - x^2"), "y")); // "x⁴ + x² − 1"

const S = ring("a, b, c, x");
String(discriminant(S.poly("a x^2 + b x + c"), "x"));              // "b² − 4ac"

const F = fractions(R);                          // the field Q(x, y, z)
String(F.parse("1/(x + y) + 1/(x - y)"));        // "2x/(x² − y²)"
String(toElementary(powerSum(R, 3)).polynomial); // "e1³ − 3e1·e2 + 3e3"
```

## What is in it

- **Rings** `ring("x, y", { field, order })` over the rationals or any RatMath
  field object (`primeField(p)`, `galoisField(p, k)`, a `NumberField`), with
  the monomial orders `lex`, `grlex`, `grevlex` or a comparison function.
- **Polynomials**: exact arithmetic, text in and out (`"x²y − (3/2)xz + 1"`,
  plain text that reads back, LaTeX), degrees, coefficients, partial
  derivatives, substitution, evaluation, change of variables.
- **Division** by one polynomial or a list (`divmod`, `rem`, `div`).
- **gcd, lcm, content, primitive part, square-free decomposition** over any
  exact field.
- **Factoring** over Q (`factor`, `isIrreducible`).
- **det, sylvester, resultant, discriminant, jacobian, hessian**.
- **Rational functions**: `fractions(R)` is the field of fractions, kept in
  lowest terms. It is a field object, so `@ratmath/linalg`'s `field.rref`,
  `field.det` and `field.inverse` work on matrices with parameters.
- **Symmetric polynomials**: `elementary`, `powerSum`, `completeHomogeneous`,
  `isSymmetric`, `toElementary`.

Coefficients handed back are Core `Rational`s over Q and the field's own
elements otherwise. Nothing is rounded.

## Limits

- Factoring uses Kronecker substitution, which is exact and simple but grows
  quickly with degree and number of variables; it throws rather than run away.
  It works over Q only.
- `squareFree` needs characteristic 0.
- Division by several polynomials gives a remainder that can depend on their
  order. Gröbner bases, which fix that, are planned (see the
  [roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md)),
  as are groups, rings and ideals.

## Install

```sh
bun add github:ratmath/structures
```

It depends on `@ratmath/core` and `@ratmath/algebraic`. Until algebraic is on
npm this repository carries it as a submodule (`vendor/algebraic`), so clone
with `--recurse-submodules`.

## Development

```sh
bun install
bun run test          # the tests in tests/
bun run docs:build    # build docs/ from docs/src/, run every example on the pages
```

Edit `docs/src/*.html` and `docs/assets/{site.js,lang.js,style.css}`; the
`docs/*.html` pages and `docs/assets/structures.js` are build output.
`docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
RatMath docs sites and should be kept in step with them.

MIT licensed.
