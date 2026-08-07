#!/usr/bin/env python3
"""Fast exact model: atom values, core values, equation satisfaction."""
import pickle, os, json, collections

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


class Model:
    def __init__(self):
        D = pickle.load(open(os.path.join(HERE, 'atoms.pkl'), 'rb'))
        self.src = D['atom_src']
        self.avars = D['atom_vars']
        self.eq_terms = D['eq_terms']
        self.eq_outer = D['eq_outer']
        self.polys = pickle.load(open(os.path.join(HERE, 'polys.pkl'), 'rb'))
        self.na = len(self.polys)
        self.ne = len(self.eq_terms)
        # atoms -> equations index
        self.atom_eqs = collections.defaultdict(list)
        for e, ts in enumerate(self.eq_terms):
            for c, a in ts:
                self.atom_eqs[a].append((e, c))
        # compiled evaluation for atoms
        self.compiled = []
        for p in self.polys:
            self.compiled.append(tuple((m, c) for m, c in p.items()))

    def atom_val(self, i, v):
        s = 0
        for m, c in self.compiled[i]:
            t = c
            for x in m:
                t *= v[x]
            s += t
        return s

    def all_atoms(self, v):
        return [self.atom_val(i, v) for i in range(self.na)]

    def core_vals(self, av):
        out = []
        for ts in self.eq_terms:
            s = 0
            for c, a in ts:
                s += c * av[a]
            out.append(s)
        return out

    def eq_fail(self, v):
        av = self.all_atoms(v)
        cv = self.core_vals(av)
        fails = []
        for e, c in enumerate(cv):
            if c == 0:
                continue
            # outer wrapper: mul by ints, squares, cubes, linear combo coefficient
            val = c
            bad = False
            for kind, k in reversed(self.eq_outer[e]):
                if kind == 'mul':
                    val = val * k
                elif kind == 'sq':
                    val = val * val
                elif kind == 'cube':
                    val = val ** 3
                elif kind == 'lin':
                    val = val * k
            if val != 0:
                fails.append(e)
        return fails, av, cv


def load_assign(path):
    d = json.load(open(path))
    v = [0] * NV
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v


if __name__ == '__main__':
    import sys, time
    M = Model()
    t = time.time()
    if len(sys.argv) > 1:
        v = load_assign(sys.argv[1])
    else:
        v = [0] * NV
    fails, av, cv = M.eq_fail(v)
    nz = [i for i, x in enumerate(av) if x != 0]
    print(f"atoms nonzero: {len(nz)}   equations failing: {len(fails)}  "
          f"score {M.ne - len(fails)}/{M.ne}  ({time.time()-t:.1f}s)")
    print("failing:", fails[:30])
