"""Agent B: streaming structural decomposition of EQUATIONS.txt into
(outer scalar) x (linear form over atoms).  Independent of the prior lab.
"""
import sys, time, pickle, collections
from bparse import parse_line, flatten_sum, split_coef, ast_key, poly

def polykey(p):
    return tuple(sorted(p.items()))

class Model:
    def __init__(self):
        self.atoms = []        # list of poly dicts
        self.atom_id = {}      # polykey -> id
        self.eqs = []          # list of (outer_scalar, [(coef, atom_id), ...], kind)
        self.kinds = collections.Counter()
        self.weird = []

    def atom(self, p):
        k = polykey(p)
        i = self.atom_id.get(k)
        if i is None:
            i = len(self.atoms)
            self.atom_id[k] = i
            self.atoms.append(p)
        return i

    def linform(self, ast):
        """flatten a sum-of-(coef*atom) into [(coef, atom_id)]"""
        out = []
        for t in flatten_sum(ast):
            c, rest = split_coef(t)
            if rest[0] == 'n':
                out.append((c * rest[1], self.atom({(): 1})))
                continue
            out.append((c, self.atom(poly(rest))))
        return out

    def add(self, idx, ast):
        # peel top-level structure
        if ast[0] == '*':
            num = 1
            nonnum = []
            for f in ast[1]:
                if f[0] == 'n':
                    num *= f[1]
                else:
                    nonnum.append(f)
            if len(nonnum) == 0:
                self.kinds['const'] += 1
                self.eqs.append((num, [], 'const'))
                return
            if len(nonnum) == 1:
                self.kinds['scalar*L'] += 1
                self.eqs.append((num, self.linform(nonnum[0]), 'scalar*L'))
                return
            keys = set(ast_key(f) for f in nonnum)
            if len(keys) == 1:
                self.kinds['scalar*L^%d' % len(nonnum)] += 1
                self.eqs.append((num, self.linform(nonnum[0]), 'pow%d' % len(nonnum)))
                return
            self.kinds['product-distinct'] += 1
            self.weird.append(idx)
            self.eqs.append((num, None, 'prod'))
            return
        if ast[0] == '+':
            # terms  c_i * S_i ; check all S_i identical
            parts = [split_coef(t) for t in ast[1]]
            keys = set(ast_key(r) for _, r in parts)
            if len(keys) == 1:
                tot = sum(c for c, _ in parts)
                self.kinds['sum-same'] += 1
                if tot == 0:
                    self.eqs.append((0, [], 'trivial'))
                else:
                    self.eqs.append((tot, self.linform(parts[0][1]), 'sum-same'))
                return
            # otherwise treat the whole thing as one linear form
            self.kinds['plain-sum'] += 1
            self.eqs.append((1, self.linform(ast), 'plain'))
            return
        self.kinds['other'] += 1
        self.weird.append(idx)
        self.eqs.append((1, self.linform(ast), 'other'))


def main():
    src = '/home/user/integer_solver/EQUATIONS.txt'
    M = Model()
    t0 = time.time()
    with open(src) as f:
        for i, line in enumerate(f):
            if not line.strip():
                continue
            M.add(i, parse_line(line))
            if (i + 1) % 5000 == 0:
                print("  %d  atoms=%d  %.1fs" % (i + 1, len(M.atoms), time.time() - t0), flush=True)
    print("eqs=%d atoms=%d  %.1fs" % (len(M.eqs), len(M.atoms), time.time() - t0))
    print("kinds:", dict(M.kinds))
    print("weird:", M.weird[:20], len(M.weird))
    with open('model.pkl', 'wb') as f:
        pickle.dump({'atoms': M.atoms, 'eqs': M.eqs}, f, -1)
    print("wrote model.pkl")

if __name__ == '__main__':
    main()
