#!/usr/bin/env python3
"""Cascade closure MOD P.  Over Z the handles absorb every quotient (verified by K9:
handles=0 seed gives 0 conflicts), so the binding content of the instance is mod p.
Given the selector bits, close the whole system mod p and count atoms that stay nonzero."""
import sys, os, json, time, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV
from circ2 import vars_of

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663


def compile_modp(n):
    o = n[0]
    if o == 'v': return 'v[%d]' % n[1]
    if o == 'c': return repr(n[1] % P)
    if o == 'neg': return '((-%s)%%P)' % compile_modp(n[1])
    return '((%s%s%s)%%P)' % (compile_modp(n[1]), o, compile_modp(n[2]))


class CascadeP:
    def __init__(self):
        self.E = Engine()
        E = self.E
        self.names = list(E.atoms.keys())
        self.nodes = [E.atoms[a] for a in self.names]
        self.avars = [sorted(vars_of(n)) for n in self.nodes]
        self.code = [compile(compile_modp(n), '<a>', 'eval') for n in self.nodes]
        self.aidx = {a: i for i, a in enumerate(self.names)}
        self.var2atoms = collections.defaultdict(list)
        for i, vs in enumerate(self.avars):
            for u in vs: self.var2atoms[u].append(i)
        self.eqrows = [[(k % P, self.aidx[a]) for k, a in row] for row in E.eqrows]
        self.n = len(self.nodes)
        self.g = {'P': P, '__builtins__': {}}

    def ev(self, i, v):
        self.g['v'] = v
        return eval(self.code[i], self.g)

    def close(self, seed, order, forbid=(), pin=None):
        """pin: {var: atom_index}.  A pinned var may ONLY be derived by its own atom.  This is
        the general forward-only guard: without it the closure happily solves a slot BACKWARDS
        from a downstream consumer (e.g. x608 = x34606*x12186 drives x12186 when the
        pass-through gate is on), which silently produces wrong slot values."""
        pin = pin or {}
        """forbid: atom indices that may NOT be used to derive a variable.  Needed to stop
        the closure running the target pin BACKWARDS into the tree instead of folding
        the leaves forward."""
        v = [0] * NV
        known = bytearray(NV)
        nunk = [len(vs) for vs in self.avars]
        done = bytearray(self.n)
        for i in forbid: done[i] = 1
        q = collections.deque(i for i in range(self.n) if nunk[i] == 1)
        derived = 0
        # provenance: trace[u] = atom index that derived u, or None if u was seeded.
        # deps[u] = the other variables that atom already knew.  Used to walk back from a
        # wrong value to the first non-forward derivation.
        self.trace = {}
        self.deps = {}

        def setv(u, val, src=None, dep=()):
            nonlocal derived
            self.trace[u] = src; self.deps[u] = tuple(dep)
            v[u] = val % P; known[u] = 1
            for j in self.var2atoms[u]:
                nunk[j] -= 1
                if nunk[j] == 1 and not done[j]: q.append(j)

        def prop():
            nonlocal derived
            while q:
                i = q.popleft()
                if done[i]: continue
                unk = [u for u in self.avars[i] if not known[u]]
                if len(unk) != 1: continue
                u = unk[0]
                if u in pin and pin[u] != i:
                    continue          # only u's own pin may derive it; leave the atom open
                v[u] = 0; c0 = self.ev(i, v)
                v[u] = 1; c1 = self.ev(i, v)
                v[u] = 2; c2 = self.ev(i, v)
                a = (c1 - c0) % P
                if (c2 - c0) % P != (2 * a) % P:
                    v[u] = 0; continue
                done[i] = 1
                if a == 0:
                    v[u] = 0; continue
                v[u] = 0
                setv(u, (-c0) * pow(a, P - 2, P) % P, src=i,
                     dep=[z for z in self.avars[i] if z != u])
                derived += 1
        prop()
        for u in order:
            if known[u]: continue
            setv(u, seed.get(u, 0)); prop()
        for u in range(NV):
            if not known[u]:
                setv(u, seed.get(u, 0)); prop()
        return v, derived

    def nzatoms(self, v):
        return [i for i in range(self.n) if self.ev(i, v)]

    def score(self, v):
        av = [self.ev(i, v) for i in range(self.n)]
        bad = 0
        for row in self.eqrows:
            t = 0
            for k, j in row: t += k * av[j]
            if t % P: bad += 1
        return bad, [i for i, x in enumerate(av) if x]


if __name__ == '__main__':
    t0 = time.time()
    C = CascadeP()
    print('build', time.time() - t0)
    d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    full = [0] * NV
    for k, val in d.items(): full[int(k[2:])] = int(val)
    vc = json.load(open(K + '/varclass.json'))
    handles, bools, others = vc['handles'], vc['bools'], vc['others']
    order = handles + bools + others + [u for u in range(NV) if u not in set(C.E.free)]
    seed = {u: 0 for u in handles}
    for u in bools: seed[u] = full[u]
    for u in others: seed[u] = full[u]
    t1 = time.time()
    v, der = C.close(seed, order)
    print('close %.2fs derived %d' % (time.time() - t1, der))
    t2 = time.time()
    bad, nz = C.score(v)
    print('score %.2fs: nonzero atoms mod p %d, failing eqs mod p %d' % (time.time() - t2, len(nz), bad))
    for i in nz[:20]: print('   ', C.names[i][:90])
