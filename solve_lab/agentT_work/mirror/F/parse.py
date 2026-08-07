#!/usr/bin/env python3
"""Independent parser for EQUATIONS.txt -> AST -> atom decomposition."""
import re, sys, json, pickle, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')

TOK = re.compile(r'\s*(x_\d+|-?\d+|\(|\)|\+|-|\*)')

class P:
    def __init__(self, s):
        self.t = []
        i = 0; n = len(s)
        while i < n:
            m = TOK.match(s, i)
            if not m: break
            self.t.append(m.group(1)); i = m.end()
        self.i = 0
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None
    def next(self):
        v = self.t[self.i]; self.i += 1; return v
    # expr := term (('+'|'-') term)*
    def expr(self):
        node = self.term()
        while self.peek() in ('+','-'):
            op = self.next(); r = self.term()
            node = (op, node, r)
        return node
    def term(self):
        node = self.atom()
        while self.peek() == '*':
            self.next(); r = self.atom()
            node = ('*', node, r)
        return node
    def atom(self):
        tk = self.peek()
        if tk == '(':
            self.next(); e = self.expr()
            assert self.next() == ')'
            return e
        if tk == '-':   # unary
            self.next(); return ('neg', self.atom())
        tk = self.next()
        if tk.startswith('x_'): return ('v', int(tk[2:]))
        return ('c', int(tk))

def parse_line(line):
    lhs = line.rsplit('=',1)[0]
    p = P(lhs)
    e = p.expr()
    assert p.i == len(p.t), (p.i, len(p.t))
    return e

def is_const(n):
    return n[0]=='c'

def flatten_sum(n, sign, out):
    """Flatten +/- tree into list of (sign, node)."""
    if n[0]=='+':
        flatten_sum(n[1], sign, out); flatten_sum(n[2], sign, out)
    elif n[0]=='-':
        flatten_sum(n[1], sign, out); flatten_sum(n[2], -sign, out)
    elif n[0]=='neg':
        flatten_sum(n[1], -sign, out)
    else:
        out.append((sign, n))

def node_str(n):
    o=n[0]
    if o=='v': return 'x%d'%n[1]
    if o=='c': return str(n[1])
    if o=='neg': return '(-%s)'%node_str(n[1])
    return '(%s%s%s)'%(node_str(n[1]), o, node_str(n[2]))

def strip_scalar(e):
    """Return (list of (coef,node)) representing the top-level structure: LHS = sum coef_i * S_i.
    Detect the repeated core S."""
    terms=[]
    flatten_sum(e, 1, terms)
    return terms

def const_val(n):
    if n[0]=='c': return n[1]
    if n[0]=='neg':
        v=const_val(n[1]); return None if v is None else -v
    if n[0]=='*':
        a=const_val(n[1]); b=const_val(n[2])
        return None if a is None or b is None else a*b
    if n[0]=='+':
        a=const_val(n[1]); b=const_val(n[2])
        return None if a is None or b is None else a+b
    if n[0]=='-':
        a=const_val(n[1]); b=const_val(n[2])
        return None if a is None or b is None else a-b
    return None

def core_of(e):
    """Return the 'S' whose vanishing is equivalent to the equation, plus a tag."""
    terms = strip_scalar(e)
    # each term: (sign, node). node may be c*S or S*S or S
    cores=[]
    for sg, nd in terms:
        # split multiplication chain
        facs=[]
        def mul(n):
            if n[0]=='*': mul(n[1]); mul(n[2])
            elif n[0]=='neg': facs.append(('c',-1)); mul(n[1])
            else: facs.append(n)
        mul(nd)
        nonc=[f for f in facs if const_val(f) is None]
        cs=[const_val(f) for f in facs if const_val(f) is not None]
        k=1
        for c in cs: k*=c
        cores.append((sg*k, nonc))
    return cores

if __name__=='__main__':
    t0=time.time()
    lines=open(EQ).read().splitlines()
    print(len(lines), 'lines')
    shapes={}
    out=[]
    for idx,ln in enumerate(lines):
        ln=ln.strip()
        if not ln: continue
        e=parse_line(ln)
        cs=core_of(e)
        out.append(cs)
        sig = tuple(len(nc) for _,nc in cs)
        shapes[sig]=shapes.get(sig,0)+1
        if idx%5000==0: print(idx, time.time()-t0, file=sys.stderr)
    print('shape signatures:', sorted(shapes.items(), key=lambda kv:-kv[1])[:20])
    pickle.dump(out, open(os.path.join(HERE,'cores.pkl'),'wb'))
    print('done', time.time()-t0)
