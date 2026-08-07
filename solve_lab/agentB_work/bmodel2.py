"""Agent B model v2: factored atoms.  An atom is a PRODUCT of polynomial factors,
so zeroing it only needs ONE factor to vanish.  Also keeps per-equation atom lists.
"""
import sys, time, pickle, collections
from bparse import parse_line, flatten_sum, split_coef, ast_key, poly

def polykey(p):
    return tuple(sorted(p.items()))

class Model:
    def __init__(self):
        self.facs = []       # list of poly dicts (factors)
        self.fac_id = {}
        self.atoms = []      # list of tuple(factor ids) sorted
        self.atom_id = {}
        self.eqs = []        # (outer, [(coef, atom_id)], kind)

    def fac(self, p):
        k = polykey(p)
        i = self.fac_id.get(k)
        if i is None:
            i = len(self.facs); self.fac_id[k] = i; self.facs.append(p)
        return i

    def atom(self, fl):
        k = tuple(sorted(fl))
        i = self.atom_id.get(k)
        if i is None:
            i = len(self.atoms); self.atom_id[k] = i; self.atoms.append(k)
        return i

    def mkatom(self, ast):
        """ast (no leading numeric coef) -> (extra_coef, atom_id)"""
        c = 1
        if ast[0] == '*':
            fl = []
            for f in ast[1]:
                if f[0] == 'n':
                    c *= f[1]
                else:
                    fl.append(self.fac(poly(f)))
            if not fl:
                return c, self.atom((self.fac({(): 1}),))
            return c, self.atom(tuple(fl))
        return c, self.atom((self.fac(poly(ast)),))

    def linform(self, ast):
        out = []
        for t in flatten_sum(ast):
            c, rest = split_coef(t)
            c2, aid = self.mkatom(rest)
            out.append((c * c2, aid))
        return out

    def add(self, ast):
        if ast[0] == '*':
            num = 1; nonnum = []
            for f in ast[1]:
                if f[0] == 'n': num *= f[1]
                else: nonnum.append(f)
            if len(nonnum) == 1:
                self.eqs.append((num, self.linform(nonnum[0]), 'sL')); return
            keys = set(ast_key(f) for f in nonnum)
            assert len(keys) == 1, "distinct product factors at top level"
            self.eqs.append((num, self.linform(nonnum[0]), 'pow%d' % len(nonnum))); return
        if ast[0] == '+':
            parts = [split_coef(t) for t in ast[1]]
            keys = set(ast_key(r) for _, r in parts)
            if len(keys) == 1:
                tot = sum(c for c, _ in parts)
                if tot == 0:
                    self.eqs.append((0, [], 'trivial')); return
                self.eqs.append((tot, self.linform(parts[0][1]), 'same')); return
            self.eqs.append((1, self.linform(ast), 'plain')); return
        self.eqs.append((1, self.linform(ast), 'other'))


def main():
    M = Model(); t0 = time.time()
    for i, line in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
        if not line.strip(): continue
        M.add(parse_line(line))
        if (i+1) % 10000 == 0:
            print("  %d facs=%d atoms=%d %.1fs" % (i+1, len(M.facs), len(M.atoms), time.time()-t0), flush=True)
    print("eqs=%d facs=%d atoms=%d %.1fs" % (len(M.eqs), len(M.facs), len(M.atoms), time.time()-t0))
    print(collections.Counter(e[2] for e in M.eqs))
    pickle.dump({'facs': M.facs, 'atoms': M.atoms, 'eqs': M.eqs}, open('model2.pkl','wb'), -1)
    print('wrote model2.pkl')

if __name__ == '__main__':
    main()
