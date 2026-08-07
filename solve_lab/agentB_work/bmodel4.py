"""Agent B model v4: gate-faithful decomposition.

Packing sums = '+' chains (only).  Each term = (int coef) * GATE.
GATE = product of factors (after peeling literals).  Factor = a polynomial
(the gate body, kept intact -- '-' and inner '+' are NOT split).
"""
import pickle, time, collections, sys
from bparse2 import parse_line, poly, ast_key

def flat_plus(a):
    if a[0] != '+': return [a]
    out = []
    for c in a[1]: out.extend(flat_plus(c))
    return out

def split_coef(a):
    if a[0] != '*': return 1, a
    c = 1; rest = []
    for f in a[1]:
        if f[0] == 'n': c *= f[1]
        else: rest.append(f)
    if not rest: return c, ('n', 1)
    if len(rest) == 1: return c, rest[0]
    return c, ('*', rest)

def polykey(p): return tuple(sorted(p.items()))

class Model:
    def __init__(self):
        self.facs = []; self.fac_id = {}
        self.atoms = []; self.atom_id = {}
        self.eqs = []
    def fac(self, p):
        k = polykey(p); i = self.fac_id.get(k)
        if i is None:
            i = len(self.facs); self.fac_id[k] = i; self.facs.append(p)
        return i
    def atom(self, fl):
        k = tuple(sorted(fl)); i = self.atom_id.get(k)
        if i is None:
            i = len(self.atoms); self.atom_id[k] = i; self.atoms.append(k)
        return i
    def mkatom(self, a):
        c = 1
        if a[0] == '*':
            fl = []
            for f in a[1]:
                if f[0] == 'n': c *= f[1]
                else: fl.append(self.fac(poly(f)))
            if not fl: return c, self.atom((self.fac({(): 1}),))
            return c, self.atom(tuple(fl))
        return c, self.atom((self.fac(poly(a)),))
    def linform(self, a):
        out = []
        for t in flat_plus(a):
            c, rest = split_coef(t)
            c2, aid = self.mkatom(rest)
            out.append((c * c2, aid))
        return out
    def add(self, a):
        if a[0] == '*':
            num = 1; nn = []
            for f in a[1]:
                if f[0] == 'n': num *= f[1]
                else: nn.append(f)
            if len(nn) == 1:
                self.eqs.append((num, self.linform(nn[0]), 'sL')); return
            ks = set(ast_key(f) for f in nn)
            if len(ks) == 1:
                self.eqs.append((num, self.linform(nn[0]), 'pow%d' % len(nn))); return
            self.eqs.append((num, None, 'PRODUCT-DISTINCT:%d' % len(nn))); return
        if a[0] == '+':
            parts = [split_coef(t) for t in a[1]]
            ks = set(ast_key(r) for _, r in parts)
            if len(ks) == 1:
                tot = sum(c for c, _ in parts)
                if tot == 0: self.eqs.append((0, [], 'trivial')); return
                self.eqs.append((tot, self.linform(parts[0][1]), 'same')); return
            self.eqs.append((1, self.linform(a), 'plain')); return
        if a[0] == '-':
            self.eqs.append((1, [(1, self.mkatom(a)[1])], 'minus')); return
        self.eqs.append((1, self.linform(a), 'other'))

def main():
    M = Model(); t0 = time.time()
    for i, line in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
        if not line.strip(): continue
        M.add(parse_line(line))
        if (i+1) % 10000 == 0:
            print("  %d facs=%d atoms=%d %.1fs" % (i+1, len(M.facs), len(M.atoms), time.time()-t0), flush=True)
    print("eqs=%d facs=%d atoms=%d %.1fs" % (len(M.eqs), len(M.facs), len(M.atoms), time.time()-t0))
    print(collections.Counter(e[2] for e in M.eqs))
    pickle.dump({'facs': M.facs, 'atoms': M.atoms, 'eqs': M.eqs}, open('model4.pkl','wb'), -1)
    print('wrote model4.pkl')

if __name__ == '__main__':
    main()
