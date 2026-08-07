#!/usr/bin/env python3
"""Agent P: fully independent recursive-descent parser for EQUATIONS.txt.

No import of any other agent's artifacts. Builds:
  - AST per equation
  - structural decomposition: eq = scalar * L^k  (verified, not assumed)
  - L = sum_j c_j * atom_j ; atom canonicalized as a polynomial dict
"""
import re, sys, json, os
from collections import defaultdict

EQ = '/home/user/integer_solver/EQUATIONS.txt'

TOK = re.compile(r'\s*(x_\d+|\d+|[()+\-*])')

def tokenize(s):
    out = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c in ' \t':
            i += 1; continue
        if c in '()+-*':
            out.append(c); i += 1; continue
        m = TOK.match(s, i)
        if not m:
            raise SyntaxError(f"bad char {c!r} at {i}: {s[max(0,i-40):i+40]!r}")
        out.append(m.group(1)); i = m.end()
    return out

# AST nodes: ('c', int) | ('v', idx) | ('+', [kids]) | ('*', [kids]) | ('-', a, b) | ('neg', a)

class P:
    def __init__(self, toks):
        self.t = toks; self.i = 0
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None
    def eat(self, x=None):
        v = self.t[self.i]
        if x is not None and v != x:
            raise SyntaxError(f"expected {x} got {v} at {self.i}")
        self.i += 1
        return v
    def expr(self):
        node = self.term()
        while self.peek() in ('+', '-'):
            op = self.eat()
            rhs = self.term()
            node = ('+', [node, rhs]) if op == '+' else ('-', node, rhs)
        return node
    def term(self):
        node = self.unary()
        while self.peek() == '*':
            self.eat()
            node = ('*', [node, self.unary()])
        return node
    def unary(self):
        if self.peek() == '-':
            self.eat(); return ('neg', self.unary())
        if self.peek() == '+':
            self.eat(); return self.unary()
        return self.atom()
    def atom(self):
        tk = self.peek()
        if tk == '(':
            self.eat('(')
            n = self.expr()
            self.eat(')')
            return n
        self.eat()
        if tk.startswith('x_'):
            return ('v', int(tk[2:]))
        return ('c', int(tk))

def parse(s):
    p = P(tokenize(s))
    n = p.expr()
    if p.i != len(p.t):
        raise SyntaxError(f"trailing tokens at {p.i}: {p.t[p.i:p.i+8]}")
    return n

# ---------- polynomial ----------
def pmul(a, b):
    r = defaultdict(int)
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            r[tuple(sorted(m1 + m2))] += c1 * c2
    return {k: v for k, v in r.items() if v}

def padd(a, b):
    r = dict(a)
    for m, c in b.items():
        r[m] = r.get(m, 0) + c
        if r[m] == 0: del r[m]
    return r

def pneg(a):
    return {m: -c for m, c in a.items()}

def topoly(n):
    t = n[0]
    if t == 'c': return {(): n[1]} if n[1] else {}
    if t == 'v': return {(n[1],): 1}
    if t == '+':
        r = {}
        for k in n[1]: r = padd(r, topoly(k))
        return r
    if t == '-': return padd(topoly(n[1]), pneg(topoly(n[2])))
    if t == 'neg': return pneg(topoly(n[1]))
    if t == '*':
        r = {(): 1}
        for k in n[1]: r = pmul(r, topoly(k))
        return r
    raise ValueError(t)

def is_const(n):
    return n[0] == 'c' or (n[0] == 'neg' and is_const(n[1]))

def constval(n):
    return n[1] if n[0] == 'c' else -constval(n[1])

# ---------- flatten additive chain into (coeff, node) summands ----------
def flat_add(n, sign=1, out=None):
    if out is None: out = []
    if n[0] == '+':
        for k in n[1]: flat_add(k, sign, out)
    elif n[0] == '-':
        flat_add(n[1], sign, out); flat_add(n[2], -sign, out)
    elif n[0] == 'neg':
        flat_add(n[1], -sign, out)
    else:
        out.append((sign, n))
    return out

def flat_mul(n, out=None):
    if out is None: out = []
    if n[0] == '*':
        for k in n[1]: flat_mul(k, out)
    elif n[0] == 'neg':
        out.append(('NEG',)); flat_mul(n[1], out)
    else:
        out.append(n)
    return out

def key(poly):
    return tuple(sorted(poly.items()))

def main():
    lines = open(EQ).read().split('\n')
    lines = [l.strip() for l in lines if l.strip()]
    print("equations:", len(lines))
    shapes = defaultdict(int)
    atom_key_to_id = {}
    atom_polys = []
    eq_rows = []      # list of (outer_scalar, power, [(coef, atom_id), ...])
    bad = []
    for ei, line in enumerate(lines):
        lhs = line.rsplit('=', 1)[0]
        ast = parse(lhs)
        # top-level: additive chain of (sign, node); each node -> mult chain
        summands = flat_add(ast)
        # for each summand, split constants from non-constant factors
        terms = []   # (scalar, tuple_of_nonconst_factor_keys, nonconst_nodes)
        for sg, nd in summands:
            fac = flat_mul(nd)
            sc = sg
            nc = []
            for f in fac:
                if f == ('NEG',): sc = -sc
                elif is_const(f): sc *= constval(f)
                else: nc.append(f)
            terms.append((sc, nc))
        # verify all summands share the same non-constant multiset (as polynomials)
        sig = None
        total = 0
        pw = None
        ncnodes = None
        ok = True
        for sc, nc in terms:
            ks = tuple(sorted(key(topoly(f)) for f in nc))
            if sig is None:
                sig = ks; pw = len(nc); ncnodes = nc
            elif ks != sig:
                ok = False
            total += sc
        if not ok:
            bad.append(ei); shapes['MIXED'] += 1
            eq_rows.append(None)
            continue
        # all non-const factors equal to each other?
        fk = set(key(topoly(f)) for f in ncnodes)
        shapes[(pw, len(fk))] += 1
        L = ncnodes[0]
        # decompose L into atoms
        parts = flat_add(L)
        row = []
        for sg, nd in parts:
            fac = flat_mul(nd)
            sc = sg; nc = []
            for f in fac:
                if f == ('NEG',): sc = -sc
                elif is_const(f): sc *= constval(f)
                else: nc.append(f)
            if len(nc) == 0:
                # pure constant summand
                ap = {(): 1}
            elif len(nc) == 1:
                ap = topoly(nc[0])
            else:
                ap = {(): 1}
                for f in nc: ap = pmul(ap, topoly(f))
            k = key(ap)
            aid = atom_key_to_id.get(k)
            if aid is None:
                aid = len(atom_polys); atom_key_to_id[k] = aid; atom_polys.append(ap)
            row.append((sc, aid))
        eq_rows.append((total, pw, len(fk), row))
    print("shapes (num_nonconst_factors, distinct_factor_polys) -> count:")
    for k, v in sorted(shapes.items(), key=lambda z: -z[1]): print("  ", k, v)
    print("mixed/bad:", len(bad), bad[:10])
    print("distinct atoms:", len(atom_polys))
    import pickle
    with open('/home/user/integer_solver/solve_lab/agentP_work/model.pkl', 'wb') as f:
        pickle.dump({'eq_rows': eq_rows, 'atom_polys': atom_polys}, f)
    print("saved model.pkl")

if __name__ == '__main__':
    main()
