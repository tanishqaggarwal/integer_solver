"""Agent B parser v3: preserves parenthesis groups ('g', child).

Discovery: the generator wraps the HEAD gate of every packing chain in an EXTRA
paren group.  So a packing chain node is ('+'|'-', [L, R]) with L = ('g', E) and E
itself a '+'/'-' node; the head gate is ('g', ('g', ...)) and stops the descent.
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
        terms = [self.term()]; ops = []
        while self.peek() in ('+','-'):
            ops.append(self.t[self.i][0]); self.i += 1
            terms.append(self.term())
        if not ops: return terms[0]
        cur = terms[0]
        for o, t in zip(ops, terms[1:]):
            cur = (o, [cur, t])
        return cur
    def term(self):
        f = [self.factor()]
        while self.peek() == '*':
            self.i += 1; f.append(self.factor())
        if len(f) == 1: return f[0]
        return ('*', f)
    def factor(self):
        k = self.peek()
        if k == '(':
            self.i += 1
            e = self.expr()
            if self.peek() != ')': raise ValueError("expected )")
            self.i += 1
            return ('g', e)
        if k == 'n':
            v = self.t[self.i][1]; self.i += 1; return ('n', v)
        if k == 'v':
            v = self.t[self.i][1]; self.i += 1; return ('v', v)
        if k == '-':
            self.i += 1; return ('*', [('n', -1), self.factor()])
        raise ValueError("unexpected %r" % (self.t[self.i:self.i+3],))

def parse_line(line):
    line = line.strip(); assert line.endswith('= 0')
    p = P(tokenize(line[:-3])); e = p.expr()
    if p.i != len(p.t): raise ValueError("trailing")
    return e

def strip_g(a):
    while a[0] == 'g': a = a[1]
    return a

def poly(a):
    k = a[0]
    if k == 'g': return poly(a[1])
    if k == 'n': return {(): a[1]} if a[1] else {}
    if k == 'v': return {(a[1],): 1}
    if k in ('+','-'):
        r = dict(poly(a[1][0])); sg = 1 if k == '+' else -1
        for m, v in poly(a[1][1]).items():
            nv = r.get(m, 0) + sg*v
            if nv: r[m] = nv
            elif m in r: del r[m]
        return r
    if k == '*':
        r = {(): 1}
        for c in a[1]:
            pc = poly(c); nr = {}
            for m1, v1 in r.items():
                for m2, v2 in pc.items():
                    mm = tuple(sorted(m1+m2))
                    nv = nr.get(mm, 0) + v1*v2
                    if nv: nr[mm] = nv
                    elif mm in nr: del nr[mm]
            r = nr
            if not r: return {}
        return r
    raise ValueError(k)

def ast_key(a):
    k = a[0]
    if k in ('n','v'): return a
    if k == 'g': return ast_key(a[1])
    return (k, tuple(ast_key(c) for c in a[1]))

def flat_pack(a):
    """Split a packing group into signed terms.  Returns list of (sign, term_ast).

    A chain link is  E = (op, [('g', E'), R])  with E' itself a '+'/'-' node.
    The HEAD gate is double-wrapped -- ('g', ('g', ...)) -- so descent stops there.
    """
    out = []
    def rec(E, sg):
        if E[0] in ('+','-'):
            L = E[1][0]; R = E[1][1]
            sg2 = sg if E[0] == '+' else -sg
            if L[0] == 'g':
                X = L[1]
                if X[0] in ('+','-'):          # chain prefix -> descend
                    rec(X, sg); out.append((sg2, R)); return
                if X[0] == 'g':                # double-wrapped HEAD gate
                    out.append((sg, L)); out.append((sg2, R)); return
        out.append((sg, E))
    rec(a[1] if a[0] == 'g' else a, 1)
    return out
