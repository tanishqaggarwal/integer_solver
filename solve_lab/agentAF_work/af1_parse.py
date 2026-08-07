#!/usr/bin/env python3
"""agent AF, step 1: independent recursive-descent parse of EQUATIONS.txt.

Produces:  af_ast.pkl  = per-equation additive decomposition into (coef, atom-AST)
Atoms are kept as canonical s-expressions (tuples) so they can be hashed.

No other agent's code or data is imported.
"""
import sys, os, pickle, time
sys.setrecursionlimit(100000)

EQ = os.path.join(os.path.dirname(__file__), '..', '..', 'EQUATIONS.txt')

# ---------------- tokenizer ----------------
def tokenize(s):
    toks = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c in '()+-*':
            toks.append(c); i += 1
        elif c.isdigit():
            j = i
            while j < n and s[j].isdigit():
                j += 1
            toks.append(('N', int(s[i:j]))); i = j
        elif c == 'x':
            assert s[i+1] == '_'
            j = i+2
            while j < n and s[j].isdigit():
                j += 1
            toks.append(('V', int(s[i+2:j]))); i = j
        elif c == ' ':
            i += 1
        else:
            raise ValueError('bad char %r at %d' % (c, i))
    return toks

# ---------------- parser ----------------
# AST node forms (tuples, hashable):
#   ('c', int)          constant
#   ('v', int)          variable id
#   ('+', a, b)         binary add   (kept BINARY, left-deep as written)
#   ('-', a, b)         binary sub
#   ('*', a, b)         binary mul   (kept BINARY, left-deep as written)
#   ('neg', a)          unary minus

class P:
    __slots__ = ('t', 'i')
    def __init__(self, t):
        self.t = t; self.i = 0
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None
    def take(self):
        x = self.t[self.i]; self.i += 1; return x

def p_expr(p):
    node = p_term(p)
    while True:
        c = p.peek()
        if c == '+':
            p.take(); node = ('+', node, p_term(p))
        elif c == '-':
            p.take(); node = ('-', node, p_term(p))
        else:
            return node

def p_term(p):
    node = p_factor(p)
    while p.peek() == '*':
        p.take(); node = ('*', node, p_factor(p))
    return node

def p_factor(p):
    c = p.peek()
    if c == '-':
        p.take(); return ('neg', p_factor(p))
    if c == '+':
        p.take(); return p_factor(p)
    if c == '(':
        p.take(); e = p_expr(p)
        assert p.take() == ')', 'unbalanced'
        return e
    tok = p.take()
    if tok[0] == 'N':
        return ('c', tok[1])
    if tok[0] == 'V':
        return ('v', tok[1])
    raise ValueError('unexpected %r' % (tok,))

def parse(s):
    p = P(tokenize(s))
    e = p_expr(p)
    assert p.i == len(p.t), 'trailing tokens'
    return e

# --------------- additive peel ---------------
def is_const(n):
    if n[0] == 'c':
        return n[1]
    if n[0] == 'neg':
        v = is_const(n[1])
        return None if v is None else -v
    return None

def chain_step(rhs):
    """A written chain step is always `(const)*(ATOM)` (or a bare constant).
    Anything else means we have reached the LEADING summand, which may itself be
    a subtraction -- peeling into it is the over-decomposition bug (agent U, §7).
    Returns (coef, atom) or None."""
    k = is_const(rhs)
    if k is not None:
        return (k, ('c', 1))
    if rhs[0] == '*':
        k = is_const(rhs[1])
        if k is not None:
            return (k, rhs[2])
        k = is_const(rhs[2])
        if k is not None:
            return (k, rhs[1])
    if rhs[0] == 'neg':
        s = chain_step(rhs[1])
        if s is not None:
            return (-s[0], s[1])
    return None

def peel_sum(node):
    """Return list of (coef:int, atomAST) reading the left-deep '+' chain.
    Peeling STOPS at the leading summand: only `(const)*(ATOM)` right operands
    are chain steps."""
    out = []
    cur = node
    while cur[0] in ('+', '-'):
        sign = 1 if cur[0] == '+' else -1
        st = chain_step(cur[2])
        if st is None:
            break
        out.append((sign * st[0], st[1]))
        cur = cur[1]
    k = is_const(cur)
    if k is not None:
        out.append((k, ('c', 1)))
    elif cur[0] == '*' and is_const(cur[1]) is not None:
        out.append((is_const(cur[1]), cur[2]))
    else:
        out.append((1, cur))
    out.reverse()
    return out

def strip_outer(node):
    """Remove a leading integer multiplier / negation wrapper: c*(BODY), (-1)*(BODY)."""
    facs = []
    cur = node
    while True:
        if cur[0] == 'neg':
            facs.append(-1); cur = cur[1]
        elif cur[0] == '*' and is_const(cur[1]) is not None:
            facs.append(is_const(cur[1])); cur = cur[2]
        elif cur[0] == '*' and is_const(cur[2]) is not None:
            facs.append(is_const(cur[2])); cur = cur[1]
        else:
            return facs, cur

def main():
    t0 = time.time()
    eqs = []
    with open(EQ) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            eqs.append(parse(line.rsplit('=', 1)[0]))
    print('parsed %d equations in %.1fs' % (len(eqs), time.time() - t0))

    # classify equation top-level shape
    from collections import Counter
    shape = Counter()
    bodies = []   # per equation: list of distinct body ASTs found at the top level
    for e in eqs:
        facs, core = strip_outer(e)
        # core may be BODY*BODY (squared), BODY, or c1*BODY + c2*BODY
        if core[0] == '*':
            l = strip_outer(core[1])[1]; r = strip_outer(core[2])[1]
            if l == r:
                shape['square']; shape['square'] += 1
                bodies.append([l])
                continue
        if core[0] in ('+', '-'):
            terms = peel_sum(core)
            uniq = set(t[1] for t in terms)
            if len(uniq) == 1:
                shape['lincomb-1body'] += 1
                bodies.append([terms[0][1]])
                continue
        shape['plain'] += 1
        bodies.append([core])
    print('top-level shapes:', dict(shape))

    with open(os.path.join(os.path.dirname(__file__), 'af_ast.pkl'), 'wb') as f:
        pickle.dump({'eqs': eqs, 'bodies': bodies}, f, 2)
    print('wrote af_ast.pkl  %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
