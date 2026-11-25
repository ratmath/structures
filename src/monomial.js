// Monomials are arrays of non-negative integer exponents, one per variable of
// the ring, and a monomial order says which of two monomials is larger.
//
// An order must be a total order in which 1 is least and which respects
// multiplication (a > b implies a·c > b·c). The three standard ones:
//
//   lex      compare exponents left to right
//   grlex    total degree first, ties broken by lex
//   grevlex  total degree first, ties broken by the LAST variable: the one
//            with the smaller last differing exponent is larger

const total = (e) => {
  let s = 0;
  for (let i = 0; i < e.length; i += 1) s += e[i];
  return s;
};

function lex(a, b) {
  for (let i = 0; i < a.length; i += 1) {
    if (a[i] !== b[i]) return a[i] - b[i];
  }
  return 0;
}

function grlex(a, b) {
  const d = total(a) - total(b);
  return d !== 0 ? d : lex(a, b);
}

function grevlex(a, b) {
  const d = total(a) - total(b);
  if (d !== 0) return d;
  for (let i = a.length - 1; i >= 0; i -= 1) {
    if (a[i] !== b[i]) return b[i] - a[i];
  }
  return 0;
}

export const ORDERS = Object.freeze({ lex, grlex, grevlex });

/** The comparison function for an order given by name or as a function (a, b) -> number. */
export function orderOf(order) {
  if (typeof order === "function") return order;
  const f = ORDERS[order];
  if (!f) throw new RangeError(`unknown monomial order "${order}"; use "lex", "grlex", "grevlex" or a comparison function`);
  return f;
}

export const totalDegree = total;
export const keyOf = (e) => e.join(",");
export const mulExp = (a, b) => a.map((x, i) => x + b[i]);
export const dividesExp = (a, b) => {
  for (let i = 0; i < a.length; i += 1) if (a[i] > b[i]) return false;
  return true;
};
export const divExp = (a, b) => a.map((x, i) => x - b[i]);
export const lcmExp = (a, b) => a.map((x, i) => Math.max(x, b[i]));
export const gcdExp = (a, b) => a.map((x, i) => Math.min(x, b[i]));

const SUP = { 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
export const superscript = (n) => String(n).replace(/\d/g, (d) => SUP[d]);
