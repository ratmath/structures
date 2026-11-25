// Reading a polynomial from text.
//
//   x^2*y - 3/2 x z + (x + y)^3
//
// + - * ^ and parentheses, with multiplication also written by juxtaposition
// ("3x", "2 x y", "(x+1)(x-1)", and "xy" when x and y are one-letter
// variables of the ring). Numbers are exact: integers, decimals such as 1.25
// (read as 5/4) and Core's repeating form 0.#3. Division is allowed by a
// constant, or when it comes out exactly. Unary minus binds looser than ^,
// so -x^2 is −(x²).

function tokenize(text) {
  const tokens = [];
  const re = /\s*(?:(\d+(?:\.\d*)?(?:#\d+)?|\.\d+(?:#\d+)?|\.#\d+)|([A-Za-z][A-Za-z0-9_]*)|(\*\*|[-+*/^()−·×])|(\S))/gy;
  let m;
  let last = 0;
  while (last < text.length && (m = re.exec(text))) {
    if (m[0].trim() === "") break;
    const at = m.index + m[0].length - (m[1] ?? m[2] ?? m[3] ?? m[4]).length;
    if (m[4] !== undefined) throw new SyntaxError(`unexpected "${m[4]}" at position ${at + 1} in "${text}"`);
    if (m[1] !== undefined) tokens.push({ type: "num", text: m[1], at });
    else if (m[2] !== undefined) tokens.push({ type: "name", text: m[2], at });
    else {
      const op = m[3] === "**" ? "^" : m[3] === "−" ? "-" : m[3] === "·" || m[3] === "×" ? "*" : m[3];
      tokens.push({ type: "op", text: op, at });
    }
    last = re.lastIndex;
  }
  return tokens;
}

// "xy" is x·y when every piece is a variable of the ring; longest names first
function splitName(vars, name) {
  if (vars.includes(name)) return [name];
  const names = [...vars].sort((a, b) => b.length - a.length);
  const parts = [];
  let rest = name;
  while (rest.length) {
    const hit = names.find((v) => rest.startsWith(v));
    if (!hit) return null;
    parts.push(hit);
    rest = rest.slice(hit.length);
  }
  return parts;
}

/** Read a polynomial of the ring R. */
export function parsePolynomial(R, text) {
  return parseWith(
    {
      vars: R.vars,
      number: (t) => R.constant(R.K.from(t)),
      variable: (v) => R.gen(v),
      // the value of a constant used as an exponent, or null
      integer: (f) => {
        if (!f.isConstant()) return null;
        const c = f.isZero() ? "0" : R.K.format(f._lead().c);
        return /^\d+$/.test(c) ? Number(c) : null;
      },
      divide: (a, b) => a.div(b),
      inexact: "(use fractions(ring) for rational functions)",
      negativePowers: false,
      name: String(R),
    },
    text,
  );
}

/**
 * The shared reader. `alg` says how to make numbers and variables; values
 * need add, sub, mul, neg, pow and isZero.
 */
export function parseWith(alg, text) {
  text = String(text);
  const tokens = tokenize(text);
  let pos = 0;
  const peek = () => tokens[pos];
  const isOp = (t, op) => t && t.type === "op" && t.text === op;
  const fail = (msg, t = peek()) => {
    throw new SyntaxError(`${msg}${t ? ` at position ${t.at + 1}` : " at the end"} in "${text}"`);
  };

  function expression() {
    let left = term();
    for (;;) {
      const t = peek();
      if (isOp(t, "+")) {
        pos += 1;
        left = left.add(term());
      } else if (isOp(t, "-")) {
        pos += 1;
        left = left.sub(term());
      } else return left;
    }
  }

  // a product; a leading sign applies to the whole product's first factor chain: -x^2 = -(x^2)
  function term() {
    let sign = 1;
    while (isOp(peek(), "-") || isOp(peek(), "+")) {
      if (peek().text === "-") sign = -sign;
      pos += 1;
    }
    let left = power();
    for (;;) {
      const t = peek();
      if (isOp(t, "*")) {
        pos += 1;
        left = left.mul(signedPower());
      } else if (isOp(t, "/")) {
        const slash = t;
        pos += 1;
        const d = signedPower();
        if (d.isZero()) fail("division by zero", slash);
        try {
          left = alg.divide(left, d);
        } catch {
          fail(`${d} does not divide ${left} exactly ${alg.inexact}`, slash);
        }
      } else if (t && (t.type === "num" || t.type === "name" || isOp(t, "("))) {
        left = left.mul(power());
      } else break;
    }
    return sign < 0 ? left.neg() : left;
  }

  function signedPower() {
    let sign = 1;
    while (isOp(peek(), "-") || isOp(peek(), "+")) {
      if (peek().text === "-") sign = -sign;
      pos += 1;
    }
    const p = power();
    return sign < 0 ? p.neg() : p;
  }

  function power() {
    const base = atom();
    if (!isOp(peek(), "^")) return base;
    const hat = peek();
    pos += 1;
    let exp;
    let negative = false;
    if (alg.negativePowers && isOp(peek(), "-")) {
      negative = true;
      pos += 1;
    }
    const t = peek();
    const whole = alg.negativePowers ? "an exponent must be a whole number" : "an exponent must be a non-negative whole number";
    if (t && t.type === "num" && /^\d+$/.test(t.text)) {
      pos += 1;
      exp = Number(t.text);
    } else if (isOp(t, "(")) {
      exp = alg.integer(atom());
      if (exp === null) fail(whole, hat);
    } else fail(whole, hat);
    if (negative) exp = -exp;
    if (isOp(peek(), "^")) fail("write a^(b^c) or (a^b)^c with parentheses");
    return base.pow(exp);
  }

  function atom() {
    const t = peek();
    if (!t) fail("the expression ends too soon");
    if (t.type === "num") {
      pos += 1;
      return alg.number(t.text);
    }
    if (t.type === "name") {
      pos += 1;
      const parts = splitName(alg.vars, t.text);
      if (!parts) fail(`"${t.text}" is not a variable of ${alg.name}`, t);
      // only the last piece takes a following exponent: xy^2 is x·y²
      let out = null;
      for (let i = 0; i < parts.length - 1; i += 1) out = out ? out.mul(alg.variable(parts[i])) : alg.variable(parts[i]);
      const lastGen = alg.variable(parts[parts.length - 1]);
      if (!out) return lastGen;
      if (isOp(peek(), "^")) {
        // put the last variable back as its own atom so power() applies the exponent to it alone
        tokens.splice(pos, 0, { type: "op", text: "*", at: t.at }, { type: "name", text: parts[parts.length - 1], at: t.at });
        return out;
      }
      return out.mul(lastGen);
    }
    if (isOp(t, "(")) {
      pos += 1;
      const inner = expression();
      if (!isOp(peek(), ")")) fail('expected ")"');
      pos += 1;
      return inner;
    }
    return fail(`unexpected "${t.text}"`);
  }

  if (!tokens.length) throw new SyntaxError("nothing to read");
  const result = expression();
  if (pos < tokens.length) fail(`unexpected "${peek().text}"`);
  return result;
}
