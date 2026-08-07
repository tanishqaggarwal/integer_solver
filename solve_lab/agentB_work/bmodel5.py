"""Agent B model v5: exact gate decomposition using paren-group boundaries."""
import pickle, time, collections, sys
from bparse3 import parse_line, poly, ast_key, strip_g, flat_pack, flat_top

def flat_prod(a):
    """Recursively flatten a product (through 'g' wrappers).
    Returns (int_coef, [non-numeric factor nodes])."""
    c = 1; fs = []
    stack = [a]
    while stack:
        n = stack.pop()
        n2 = strip_g(n)
        if n2[0] == 'n': c *= n2[1]
        elif n2[0] == '*': stack.extend(n2[1])
        else: fs.append(n)
    return c, fs

def split_coef(a):
    a = strip_g(a)
    if a[0] != '*': return 1, a
    c = 1; rest = []
    for f in a[1]:
        f2 = strip_g(f)
        if f2[0] == 'n': c *= f2[1]
        else: rest.append(f)
    if not rest: return c, ('n', 1)
    if len(rest) == 1: return c, rest[0]
    return c, ('*', rest)

def polykey(p): return tuple(sorted(p.items()))

class Model:
    def __init__(self):
        self.facs = []; self.fac_id = {}
        self.atoms = []; self.atom_id = {}
        self.eqs = []; self.bad = []
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
        c, fs = flat_prod(a)
        if not fs: return c, self.atom((self.fac({(): 1}),))
        return c, self.atom(tuple(self.fac(poly(f)) for f in fs))
    def linform(self, a):
        out = []
        for sg, t in flat_pack(a):
            c2, aid = self.mkatom(t)
            out.append((sg * c2, aid))
        return out
    def add(self, i, a):
        a2 = strip_g(a)
        if a2[0] not in ('+','-'):
            num, nn = flat_prod(a)
            if len(nn) >= 1:
                ks = set(ast_key(f) for f in nn)
                if len(ks) == 1:
                    self.eqs.append((num, self.linform(nn[0]),
                                     'sL' if len(nn) == 1 else 'pow%d' % len(nn))); return
                self.bad.append(i); self.eqs.append((num, None, 'PD%d' % len(nn))); return
            self.eqs.append((num, [], 'const')); return
        if True:
            parts = []
            for sg, t in flat_top(a2):
                c, fs = flat_prod(t)
                parts.append((sg*c, tuple(sorted(ast_key(f) for f in fs)), fs))
            ks = set(k for _, k, _ in parts)
            if len(ks) == 1:
                tot = sum(c for c, _, _ in parts)
                if tot == 0: self.eqs.append((0, [], 'trivial')); return
                fs = parts[0][2]
                if len(set(ast_key(f) for f in fs)) == 1:
                    self.eqs.append((tot, self.linform(fs[0]),
                                     'same' if len(fs) == 1 else 'same_pow%d' % len(fs))); return
                self.bad.append(i); self.eqs.append((tot, None, 'samePD')); return
            self.eqs.append((1, self.linform(a2), 'plain')); return
        self.bad.append(i); self.eqs.append((1, self.linform(a2), 'other'))

def main():
    M = Model(); t0 = time.time()
    for i, line in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
        if not line.strip(): continue
        M.add(i, parse_line(line))
        if (i+1) % 10000 == 0:
            print("  %d facs=%d atoms=%d %.1fs" % (i+1, len(M.facs), len(M.atoms), time.time()-t0), flush=True)
    print("eqs=%d facs=%d atoms=%d %.1fs" % (len(M.eqs), len(M.facs), len(M.atoms), time.time()-t0))
    print(collections.Counter(e[2] for e in M.eqs)); print('bad', M.bad[:10], len(M.bad))
    pickle.dump({'facs': M.facs, 'atoms': M.atoms, 'eqs': M.eqs}, open('model5.pkl','wb'), -1)
    print('wrote model5.pkl')

if __name__ == '__main__':
    main()
