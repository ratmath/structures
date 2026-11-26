# Changelog

## Unreleased

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
