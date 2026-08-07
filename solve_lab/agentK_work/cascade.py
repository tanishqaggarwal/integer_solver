#!/usr/bin/env python3
"""Generic exact-integer cascade closer over ALL 39,033 atoms.

No DAG orientation is assumed. An atom is usable when exactly one of its variables is
still unknown; we then test linearity in that variable by 3-point evaluation and solve
over Z (only if the coefficient divides exactly).  Everything is exact integer.
"""
import sys, os, json, time, collections, pickle
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV, compile_node
from circ2 import vars_of

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


class Cascade:
    def __init__(self, E=None):
        self.E = E or Engine()
        E = self.E
        self.atomnames = list(E.atoms.keys())
        self.nodes = [E.atoms[a] for a in self.atomnames]
        self.avars = [sorted(vars_of(n)) for n in self.nodes]
        self.code = [compile(compile_node(n), '<a>', 'eval') for n in self.nodes]
        self.aidx = {a: i for i, a in enumerate(self.atomnames)}
        self.var2atoms = collections.defaultdict(list)
        for i, vs in enumerate(self.avars):
            for u in vs:
                self.var2atoms[u].append(i)
        # equation rows over all atoms
        self.eqrows = [[(k, self.aidx[a]) for k, a in row] for row in E.eqrows]

    def evala(self, i, v):
        return eval(self.code[i], {'v': v, '__builtins__': {}})

    def close(self, seed, verbose=True):
        """seed: dict var->value. Returns (v, known, conflicts)."""
        v = [0] * NV
        known = bytearray(NV)
        for k, val in seed.items():
            v[k] = val; known[k] = 1
        nunk = [sum(1 for u in vs if not known[u]) for vs in self.avars]
        queue = collections.deque(i for i in range(len(self.nodes)) if nunk[i] == 1)
        conflicts = []
        done = bytearray(len(self.nodes))
        t0 = time.time(); nassign = 0
        while queue:
            i = queue.popleft()
            if done[i]: continue
            unk = [u for u in self.avars[i] if not known[u]]
            if len(unk) != 1:
                continue
            u = unk[0]
            v[u] = 0; c0 = self.evala(i, v)
            v[u] = 1; c1 = self.evala(i, v)
            v[u] = 2; c2 = self.evala(i, v)
            a = c1 - c0
            if c2 - c0 != 2 * a:
                v[u] = 0
                continue  # nonlinear in u; leave for later
            done[i] = 1
            if a == 0:
                v[u] = 0
                if c0 != 0: conflicts.append(('unsat-const', self.atomnames[i], c0))
                continue
            if (-c0) % a != 0:
                v[u] = 0
                conflicts.append(('nondiv', self.atomnames[i], u, a, c0))
                continue
            v[u] = (-c0) // a
            known[u] = 1; nassign += 1
            for j in self.var2atoms[u]:
                nunk[j] -= 1
                if nunk[j] == 1 and not done[j]: queue.append(j)
        if verbose:
            print('cascade: assigned %d, known %d/%d, conflicts %d, %.1fs'
                  % (nassign, sum(known), NV, len(conflicts), time.time() - t0))
        return v, known, conflicts

    def residuals(self, v):
        return [i for i in range(len(self.nodes)) if self.evala(i, v)]

    def score(self, v):
        av = [self.evala(i, v) for i in range(len(self.nodes))]
        bad = []
        for e, row in enumerate(self.eqrows):
            t = 0
            for k, j in row: t += k * av[j]
            if t: bad.append(e)
        nz = [i for i, x in enumerate(av) if x]
        return bad, nz


if __name__ == '__main__':
    C = Cascade()
    d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    full = [0] * NV
    for k, val in d.items(): full[int(k[2:])] = int(val)
    freeset = set(C.E.free)
    # seed: every free input whose deliverable value is 0 or 1 (the selector layer)
    seed = {u: full[u] for u in freeset if full[u] in (0, 1)}
    print('seed size', len(seed), 'of', len(freeset), 'free inputs')
    v, known, conf = C.close(seed)
    bad, nz = C.score(v)
    print('nonzero atoms', len(nz), 'failing eqs', len(bad), '=> score', 39033 - len(bad))
    for c in conf[:20]: print('  conflict', c[0], c[1][:80])
    print('unknown vars remaining', NV - sum(known))
    json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]},
              open('/home/user/integer_solver/solve_lab/agentK_work/cand_cascade.json', 'w'))
