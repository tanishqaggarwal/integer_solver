"""Agent B: independent parser for EQUATIONS.txt.

Grammar (fully parenthesised, but we parse generally):
  expr := term (('+'|'-') term)*
  term := factor ('*' factor)*
  factor := '(' expr ')' | INT | VAR | '-' factor

AST nodes (tuples):
  ('n', int)          integer literal
  ('v', int)          variable index
  ('+', [children])   sum
  ('*', [children])   product
"""
import re, sys, pickle, time

TOK = re.compile(r'\s*(?:(x_(\d+))|(-?\d+)|([()+*-]))')

def tokenize(s):
    toks = []
    pos = 0
    n = len(s)
    ap = toks.append
    while pos < n:
        m = TOK.match(s, pos)
        if not m:
            if s[pos:].strip() == '':
                break
            raise ValueError("tok fail at %d: %r" % (pos, s[pos:pos+40]))
        pos = m.end()
        if m.group(1) is not None:
            ap(('v', int(m.group(2))))
        elif m.group(3) is not None:
            ap(('n', int(m.group(3))))
        else:
            ap((m.group(4), None))
    return toks


class P:
    __slots__ = ('t', 'i')
    def __init__(self, toks):
        self.t = toks
        self.i = 0
    def peek(self):
        return self.t[self.i][0] if self.i < len(self.t) else None
    def expr(self):
        # note: tokenizer already folds '-5' into a literal when preceded by nothing,
        # but binary minus appears as separate token in patterns like (a)-(b)
        terms = [self.term()]
        while True:
            k = self.peek()
            if k == '+':
                self.i += 1
                terms.append(self.term())
            elif k == '-':
                self.i += 1
                terms.append(('*', [('n', -1), self.term()]))
            elif k == 'n' and self.t[self.i][1] < 0:
                # e.g. "(a) (-5)" cannot happen; but "(a)-5" tokenizes as n(-5)
                terms.append(('n', self.t[self.i][1]))
                self.i += 1
            else:
                break
        if len(terms) == 1:
            return terms[0]
        return ('+', terms)
    def term(self):
        f = [self.factor()]
        while self.peek() == '*':
            self.i += 1
            f.append(self.factor())
        if len(f) == 1:
            return f[0]
        return ('*', f)
    def factor(self):
        k = self.peek()
        if k == '(':
            self.i += 1
            e = self.expr()
            if self.peek() != ')':
                raise ValueError("expected ) at %d" % self.i)
            self.i += 1
            return e
        if k == 'n':
            v = self.t[self.i][1]; self.i += 1
            return ('n', v)
        if k == 'v':
            v = self.t[self.i][1]; self.i += 1
            return ('v', v)
        if k == '-':
            self.i += 1
            return ('*', [('n', -1), self.factor()])
        raise ValueError("unexpected token %r at %d" % (self.t[self.i:self.i+3], self.i))


def parse_line(line):
    line = line.strip()
    assert line.endswith('= 0'), line[-20:]
    body = line[:-3]
    p = P(tokenize(body))
    e = p.expr()
    if p.i != len(p.t):
        raise ValueError("trailing tokens at %d/%d" % (p.i, len(p.t)))
    return e


# ---------- polynomial normal form ----------
# poly = dict: monomial (sorted tuple of var indices, with repetition) -> int coeff

def poly(ast):
    k = ast[0]
    if k == 'n':
        return {(): ast[1]} if ast[1] else {}
    if k == 'v':
        return {(ast[1],): 1}
    if k == '+':
        r = {}
        for c in ast[1]:
            for m, v in poly(c).items():
                nv = r.get(m, 0) + v
                if nv:
                    r[m] = nv
                elif m in r:
                    del r[m]
        return r
    if k == '*':
        r = {(): 1}
        for c in ast[1]:
            pc = poly(c)
            nr = {}
            for m1, v1 in r.items():
                for m2, v2 in pc.items():
                    m = tuple(sorted(m1 + m2))
                    nv = nr.get(m, 0) + v1 * v2
                    if nv:
                        nr[m] = nv
                    elif m in nr:
                        del nr[m]
            r = nr
            if not r:
                return {}
        return r
    raise ValueError(k)


def flatten_sum(ast):
    """Flatten nested '+' into a flat list of children."""
    if ast[0] != '+':
        return [ast]
    out = []
    for c in ast[1]:
        out.extend(flatten_sum(c))
    return out


def split_coef(ast):
    """Return (int_coef, rest_ast) for a product node with leading literals."""
    if ast[0] != '*':
        return 1, ast
    c = 1
    rest = []
    for f in ast[1]:
        if f[0] == 'n':
            c *= f[1]
        else:
            rest.append(f)
    if not rest:
        return c, ('n', 1)
    if len(rest) == 1:
        return c, rest[0]
    return c, ('*', rest)


def ast_key(ast):
    """Canonical hashable key for an AST (structural, no algebra)."""
    k = ast[0]
    if k in ('n', 'v'):
        return ast
    return (k, tuple(ast_key(c) for c in ast[1]))


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else '/home/user/integer_solver/EQUATIONS.txt'
    out = sys.argv[2] if len(sys.argv) > 2 else 'asts.pkl'
    t0 = time.time()
    asts = []
    with open(src) as f:
        for i, line in enumerate(f):
            if not line.strip():
                continue
            try:
                asts.append(parse_line(line))
            except Exception as ex:
                print("LINE %d FAIL: %s" % (i, ex))
                print(line[:200])
                raise
            if (i + 1) % 5000 == 0:
                print("  %d lines  %.1fs" % (i + 1, time.time() - t0), flush=True)
    print("parsed %d equations in %.1fs" % (len(asts), time.time() - t0))
    with open(out, 'wb') as f:
        pickle.dump(asts, f, -1)
    print("wrote", out)


if __name__ == '__main__':
    main()
