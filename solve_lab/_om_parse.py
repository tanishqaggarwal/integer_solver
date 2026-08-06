#!/usr/bin/env python3
"""Off-manifold: parse EQUATIONS.txt into  eq = M * (sum_j c_j * Atom_j)^k .
Leaf 'atoms' are maximal subexpressions that are not (+) combinations or
const-multiples.  Produces atom<->equation incidence with exact integer coeffs."""
import re, sys, json, pickle, os

TOK = re.compile(r'\s*(\(|\)|\*|\+|-|x_\d+|\d+)')

def tokenize(s):
    s = s.strip()
    out = []; i = 0; n = len(s)
    while i < n:
        if s[i].isspace(): i += 1; continue
        m = TOK.match(s, i)
        if not m: raise ValueError('tok fail at %d: %r' % (i, s[i:i+40]))
        out.append(m.group(1)); i = m.end()
    return out

# AST nodes: ('num',v) ('var',i) ('add',a,b) ('sub',a,b) ('mul',a,b) ('neg',a)
class P:
    def __init__(self, toks): self.t = toks; self.i = 0
    def peek(self): return self.t[self.i] if self.i < len(self.t) else None
    def eat(self, x=None):
        v = self.t[self.i]; self.i += 1
        if x is not None and v != x: raise ValueError('want %s got %s' % (x, v))
        return v
    def expr(self):
        n = self.term()
        while self.peek() in ('+', '-'):
            op = self.eat()
            r = self.term()
            n = ('add', n, r) if op == '+' else ('sub', n, r)
        return n
    def term(self):
        n = self.unary()
        while self.peek() == '*':
            self.eat('*'); n = ('mul', n, self.unary())
        return n
    def unary(self):
        if self.peek() == '-':
            self.eat('-'); return ('neg', self.unary())
        return self.atom()
    def atom(self):
        tk = self.peek()
        if tk == '(':
            self.eat('('); n = self.expr(); self.eat(')'); return n
        tk = self.eat()
        if tk.startswith('x_'): return ('var', int(tk[2:]))
        return ('num', int(tk))

def key(n):
    """canonical string of a subtree (for structural equality + atom identity)"""
    t = n[0]
    if t == 'num': return str(n[1])
    if t == 'var': return 'x%d' % n[1]
    if t == 'neg': return '-(%s)' % key(n[1])
    op = {'add': '+', 'sub': '-', 'mul': '*'}[t]
    return '(%s%s%s)' % (key(n[1]), op, key(n[2]))

def constval(n):
    """return int if subtree is a pure constant else None"""
    t = n[0]
    if t == 'num': return n[1]
    if t == 'var': return None
    if t == 'neg':
        a = constval(n[1]); return None if a is None else -a
    a = constval(n[1]); b = constval(n[2])
    if a is None or b is None: return None
    return a + b if t == 'add' else (a - b if t == 'sub' else a * b)

def is_cmul(n):
    """True if n is (const)*(something non-const) -- the body-level term shape"""
    if n[0] != 'mul': return False
    ca, cb = constval(n[1]), constval(n[2])
    return (ca is not None) != (cb is not None)

def linearize(n, out, c=1):
    """decompose into sum_j c_j * leaf_j ; out is dict key->coeff, also records leaf ast.
    A body is a LEFT-NESTED chain  ((A_0 + c1*A_1) + c2*A_2) + ...  so an 'add' node is a
    body-level split ONLY when its right child has the shape (const)*(expr).  Otherwise the
    '+' belongs inside a single atom (e.g. (x_a*x_b) + x_c)."""
    t = n[0]
    if t == 'add' and is_cmul(n[2]):
        linearize(n[1], out, c); linearize(n[2], out, c); return
    if t == 'sub':
        # only treat as body-level split if BOTH sides are non-constant sums? no:
        # atoms use '-' internally.  Treat 'sub' as a LEAF (atom) unless it is the
        # top level combination, which by construction always uses '+'.
        pass
    if t == 'neg':
        linearize(n[1], out, -c); return
    if t == 'mul':
        a, b = n[1], n[2]
        ca, cb = constval(a), constval(b)
        if ca is not None and cb is not None:
            return  # pure constant term (shouldn't happen at body level)
        if ca is not None:
            linearize(b, out, c * ca); return
        if cb is not None:
            linearize(a, out, c * cb); return
        if key(a) == key(b):
            linearize(a, out, c); return  # perfect power  -> body == 0
        # genuine monomial: leaf
    k = key(n)
    out.setdefault(k, [0, n])
    out[k][0] += c

def parse_eq(line):
    lhs = line.rsplit('=', 1)[0]
    ast = P(tokenize(lhs)).expr()
    out = {}
    linearize(ast, out)
    return {k: v[0] for k, v in out.items() if v[0] != 0}, {k: v[1] for k, v in out.items()}

def evalast(n, V):
    t = n[0]
    if t == 'num': return n[1]
    if t == 'var': return V[n[1]]
    if t == 'neg': return -evalast(n[1], V)
    a = evalast(n[1], V); b = evalast(n[2], V)
    return a + b if t == 'add' else (a - b if t == 'sub' else a * b)

def astvars(n, s=None):
    if s is None: s = set()
    if n[0] == 'var': s.add(n[1])
    elif n[0] == 'neg': astvars(n[1], s)
    elif n[0] in ('add', 'sub', 'mul'): astvars(n[1], s); astvars(n[2], s)
    return s

if __name__ == '__main__':
    lines = [L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
    print('n eq', len(lines))
    eqatoms = []   # list of dict atomkey -> coeff
    astof = {}
    for i, L in enumerate(lines):
        d, asts = parse_eq(L)
        eqatoms.append(d)
        for k, a in asts.items():
            if k not in astof: astof[k] = a
        if i % 5000 == 0: print(' ', i, file=sys.stderr)
    print('distinct leaf atoms:', len(astof))
    with open('_om_parsed2.pkl', 'wb') as f:
        pickle.dump({'eqatoms': eqatoms, 'astof': astof}, f)
    print('saved')
