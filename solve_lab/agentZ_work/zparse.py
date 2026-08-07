#!/usr/bin/env python3
"""Agent Z: independent recursive-descent parse of EQUATIONS.txt.

Goal: decompose every equation into  scalar * L^k = 0  with L a linear
combination of "atoms", then find every atom / every equation that mentions
two or more distinct SELECTOR variables.

No code from any other agent dir is read or imported.
"""
import sys, re, os, pickle, random

EQ = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'EQUATIONS.txt')

# ---------------------------------------------------------------- tokenizer
TOK = re.compile(r'\s*(\(|\)|\+|-|\*|x_\d+|\d+)')

def tokenize(s):
    out = []
    i = 0
    n = len(s)
    while i < n:
        m = TOK.match(s, i)
        if not m:
            if s[i].isspace():
                i += 1
                continue
            raise ValueError("bad char %r at %d" % (s[i], i))
        out.append(m.group(1))
        i = m.end()
    return out

# ---------------------------------------------------------------- AST
# node forms:
#   ('n', int)              numeric literal
#   ('v', int)              variable id
#   ('+', (t1,t2,...))      sum
#   ('*', (f1,f2,...))      product
#   ('neg', t)              unary minus

class P:
    def __init__(self, toks):
        self.t = toks
        self.i = 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None

    def next(self):
        v = self.t[self.i]
        self.i += 1
        return v

    def expr(self):
        terms = []
        sign = 1
        if self.peek() == '-':
            self.next(); sign = -1
        elif self.peek() == '+':
            self.next()
        t = self.term()
        terms.append(t if sign == 1 else neg(t))
        while self.peek() in ('+', '-'):
            op = self.next()
            t = self.term()
            terms.append(t if op == '+' else neg(t))
        if len(terms) == 1:
            return terms[0]
        return ('+', tuple(terms))

    def term(self):
        fs = [self.factor()]
        while self.peek() == '*':
            self.next()
            fs.append(self.factor())
        if len(fs) == 1:
            return fs[0]
        return ('*', tuple(fs))

    def factor(self):
        tk = self.peek()
        if tk == '(':
            self.next()
            e = self.expr()
            assert self.next() == ')'
            return e
        if tk == '-':
            self.next()
            return neg(self.factor())
        tk = self.next()
        if tk.startswith('x_'):
            return ('v', int(tk[2:]))
        return ('n', int(tk))

def neg(t):
    if t[0] == 'n':
        return ('n', -t[1])
    return ('*', (('n', -1), t))

def parse(s):
    p = P(tokenize(s))
    e = p.expr()
    assert p.i == len(p.t), "trailing tokens"
    return e

# ---------------------------------------------------------------- helpers
def flat_add(t):
    """flatten nested sums into a term list"""
    if t[0] != '+':
        return [t]
    out = []
    for s in t[1]:
        out.extend(flat_add(s))
    return out

def flat_mul(t):
    if t[0] != '*':
        return [t]
    out = []
    for s in t[1]:
        out.extend(flat_mul(s))
    return out

def split_coeff(t):
    """term -> (integer coeff, core node) stripping numeric factors"""
    fs = flat_mul(t)
    c = 1
    core = []
    for f in fs:
        if f[0] == 'n':
            c *= f[1]
        else:
            core.append(f)
    if not core:
        return c, None
    if len(core) == 1:
        return c, core[0]
    return c, ('*', tuple(core))

def key(t):
    """structural canonical key (order preserved -- purely syntactic)"""
    if t[0] == 'n':
        return 'n%d' % t[1]
    if t[0] == 'v':
        return 'x%d' % t[1]
    return '(' + t[0] + ''.join(key(s) for s in t[1]) + ')'

def varset(t, acc=None):
    if acc is None:
        acc = set()
    st = [t]
    while st:
        u = st.pop()
        if u[0] == 'v':
            acc.add(u[1])
        elif u[0] in '+*':
            st.extend(u[1])
    return acc

# ---------------------------------------------------------------- reduce to linear form L
def reduce_L(t, depth=0):
    """Strip scalar multiples and powers: return the node L such that the
    equation is  scalar * L^k = 0 .  Returns (L, note)."""
    if depth > 20:
        return t, 'depthcap'
    if t[0] == 'v':
        return t, 'var'
    if t[0] == 'n':
        return t, 'const'
    if t[0] == '*':
        c, core = split_coeff(t)
        if core is None:
            return t, 'const'
        fs = flat_mul(core)
        ks = set(key(f) for f in fs)
        if len(ks) == 1:
            return reduce_L(fs[0], depth + 1)
        # genuinely different factors: this is the linear form itself
        return t, 'mulmixed'
    # sum
    terms = flat_add(t)
    cores = []
    for tm in terms:
        c, core = split_coeff(tm)
        if core is None:
            cores.append(None)
        else:
            cores.append(core)
    nonnull = [c for c in cores if c is not None]
    if len(terms) <= 2 and len(nonnull) == len(terms) and len(set(key(c) for c in nonnull)) == 1:
        # C1*L + C2*L  shape
        return reduce_L(nonnull[0], depth + 1)
    return t, 'sum'

def atoms_of(L):
    """Split the linear form into (coeff, atom_node) list."""
    if L[0] != '+':
        c, core = split_coeff(L)
        return [(c, core if core is not None else ('n', 1))]
    out = []
    for tm in flat_add(L):
        c, core = split_coeff(tm)
        out.append((c, core if core is not None else ('n', 1)))
    return out

# ---------------------------------------------------------------- numeric eval (validation)
def ev(t, env):
    if t[0] == 'n':
        return t[1]
    if t[0] == 'v':
        return env.get(t[1], 0)
    if t[0] == '+':
        s = 0
        for u in t[1]:
            s += ev(u, env)
        return s
    p = 1
    for u in t[1]:
        p *= ev(u, env)
    return p

# ---------------------------------------------------------------- main
def main():
    sys.setrecursionlimit(100000)
    lines = []
    with open(EQ) as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line.rsplit('=', 1)[0])
    print("equations:", len(lines))

    eq_atoms = []     # per equation: list of (coeff, atomkey)
    atom_tab = {}     # atomkey -> node
    notes = {}
    for i, lhs in enumerate(lines):
        E = parse(lhs)
        L, note = reduce_L(E)
        notes[note] = notes.get(note, 0) + 1
        al = atoms_of(L)
        ea = []
        for c, a in al:
            k = key(a)
            if k not in atom_tab:
                atom_tab[k] = a
            ea.append((c, k))
        eq_atoms.append(ea)
        if i % 5000 == 0:
            print("  ...", i, flush=True)
    print("reduce notes:", notes)
    print("distinct atoms:", len(atom_tab))

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zmodel.pkl'), 'wb') as fh:
        pickle.dump({'eq_atoms': eq_atoms, 'atom_tab': atom_tab}, fh)
    print("saved zmodel.pkl")

if __name__ == '__main__':
    main()
