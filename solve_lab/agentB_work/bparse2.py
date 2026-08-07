"""Agent B parser v2: keeps binary '-' distinct from explicit (literal)*(...) products,
and flattens products.  This lets us recover the ORIGINAL gate grouping:
the packing sums are chains of '+' whose terms are (literal)*(gate);
gate bodies use '-' and bare '+'.

AST:
  ('n', k) ('v', i)
  ('+', [children])   from '+' operators only
  ('-', [a, b])       from a binary '-' operator
  ('*', [factors])    flattened product
"""
import re

TOK = re.compile(r'\s*(?:(x_(\d+))|(\d+)|([()+*-]))')

def tokenize(s):
    toks = []; pos = 0; n = len(s); ap = toks.append
    while pos < n:
        m = TOK.match(s, pos)
        if not m:
            if not s[pos:].strip(): break
            raise ValueError("tok fail %r" % s[pos:pos+40])
        pos = m.end()
        if m.group(1) is not None: ap(('v', int(m.group(2))))
        elif m.group(3) is not None: ap(('n', int(m.group(3))))
        else: ap((m.group(4), None))
    return toks

class P:
    __slots__ = ('t','i')
    def __init__(self, t): self.t = t; self.i = 0
    def peek(self): return self.t[self.i][0] if self.i < len(self.t) else None
    def expr(self):
        terms = [self.term()]
        ops = []
        while True:
            k = self.peek()
            if k in ('+', '-'):
                self.i += 1
                ops.append(k)
                terms.append(self.term())
            else:
                break
        if not ops:
            return terms[0]
        if all(o == '+' for o in ops):
            return ('+', terms)
        # left-associative mixed chain
        cur = terms[0]
        for o, t in zip(ops, terms[1:]):
            cur = ('-', [cur, t]) if o == '-' else ('+', [cur, t])
        return cur
    def term(self):
        f = [self.factor()]
        while self.peek() == '*':
            self.i += 1
            f.append(self.factor())
        if len(f) == 1: return f[0]
        out = []
        for g in f:
            if g[0] == '*': out.extend(g[1])
            else: out.append(g)
        return ('*', out)
    def factor(self):
        k = self.peek()
        if k == '(':
            self.i += 1
            e = self.expr()
            if self.peek() != ')': raise ValueError("expected )")
            self.i += 1
            return e
        if k == 'n':
            v = self.t[self.i][1]; self.i += 1; return ('n', v)
        if k == 'v':
            v = self.t[self.i][1]; self.i += 1; return ('v', v)
        if k == '-':
            self.i += 1
            return ('*', [('n', -1), self.factor()])
        raise ValueError("unexpected %r" % (self.t[self.i:self.i+3],))

def parse_line(line):
    line = line.strip()
    assert line.endswith('= 0')
    p = P(tokenize(line[:-3]))
    e = p.expr()
    if p.i != len(p.t): raise ValueError("trailing")
    return e

def poly(ast):
    k = ast[0]
    if k == 'n': return {(): ast[1]} if ast[1] else {}
    if k == 'v': return {(ast[1],): 1}
    if k == '-':
        a = poly(ast[1][0]); b = poly(ast[1][1])
        r = dict(a)
        for m, v in b.items():
            nv = r.get(m, 0) - v
            if nv: r[m] = nv
            elif m in r: del r[m]
        return r
    if k == '+':
        r = {}
        for c in ast[1]:
            for m, v in poly(c).items():
                nv = r.get(m, 0) + v
                if nv: r[m] = nv
                elif m in r: del r[m]
        return r
    if k == '*':
        r = {(): 1}
        for c in ast[1]:
            pc = poly(c); nr = {}
            for m1, v1 in r.items():
                for m2, v2 in pc.items():
                    m = tuple(sorted(m1+m2))
                    nv = nr.get(m, 0) + v1*v2
                    if nv: nr[m] = nv
                    elif m in nr: del nr[m]
            r = nr
            if not r: return {}
        return r
    raise ValueError(k)

def ast_key(a):
    k = a[0]
    if k in ('n','v'): return a
    return (k, tuple(ast_key(c) for c in a[1]))
