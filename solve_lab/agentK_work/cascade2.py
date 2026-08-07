#!/usr/bin/env python3
"""Incremental cascade: propagate first, only guess a variable when nothing else is forced."""
import sys, os, json, time, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P


class Inc:
    def __init__(self, C):
        self.C = C
        self.n = len(C.nodes)

    def reset(self):
        C = self.C
        self.v = [0] * NV
        self.known = bytearray(NV)
        self.nunk = [len(vs) for vs in C.avars]
        self.done = bytearray(self.n)
        self.queue = collections.deque(i for i in range(self.n) if self.nunk[i] == 1)
        self.conflicts = []
        self.nassign = 0
        self.derived = []

    def set(self, u, val):
        if self.known[u]:
            return self.v[u] == val
        self.v[u] = val; self.known[u] = 1
        for j in self.C.var2atoms[u]:
            self.nunk[j] -= 1
            if self.nunk[j] == 1 and not self.done[j]: self.queue.append(j)
        return True

    def propagate(self):
        C = self.C; v = self.v
        while self.queue:
            i = self.queue.popleft()
            if self.done[i]: continue
            unk = [u for u in C.avars[i] if not self.known[u]]
            if len(unk) != 1: continue
            u = unk[0]
            v[u] = 0; c0 = C.evala(i, v)
            v[u] = 1; c1 = C.evala(i, v)
            v[u] = 2; c2 = C.evala(i, v)
            a = c1 - c0
            if c2 - c0 != 2 * a:
                v[u] = 0
                continue
            self.done[i] = 1
            if a == 0:
                v[u] = 0
                if c0 != 0: self.conflicts.append(('unsat', C.atomnames[i], c0))
                continue
            if (-c0) % a != 0:
                v[u] = 0
                self.conflicts.append(('nondiv', C.atomnames[i], u, a, c0))
                continue
            val = (-c0) // a
            v[u] = 0
            self.set(u, val)
            self.nassign += 1
            self.derived.append(u)

    def run(self, seedvals, order, verbose=True):
        """seedvals: dict var->value used only when a var is still unknown and unforced."""
        self.reset()
        self.propagate()
        t0 = time.time()
        for u in order:
            if self.known[u]: continue
            self.set(u, seedvals.get(u, 0))
            self.propagate()
        # anything still unknown -> 0
        rest = [u for u in range(NV) if not self.known[u]]
        for u in rest:
            if not self.known[u]:
                self.set(u, seedvals.get(u, 0)); self.propagate()
        if verbose:
            print('inc cascade: derived %d, conflicts %d, %.1fs' % (self.nassign, len(self.conflicts), time.time() - t0))
        return self.v


if __name__ == '__main__':
    C = Cascade()
    d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    full = [0] * NV
    for k, val in d.items(): full[int(k[2:])] = int(val)
    I = Inc(C)
    order = sorted(range(NV))
    v = I.run({u: full[u] for u in range(NV)}, order)
    bad, nz = C.score(v)
    print('nonzero atoms', len(nz), 'failing eqs', len(bad), '=> score', 39033 - len(bad))
    for c in I.conflicts[:20]: print('  conflict', c[0], str(c[1])[:80])
    for i in nz[:20]: print('  nz atom', C.atomnames[i][:90])
    json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]},
              open(K + '/cand_cascade2.json', 'w'))
    print('n differing from deliverable:', sum(1 for i in range(NV) if v[i] != full[i]))
