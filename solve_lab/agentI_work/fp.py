#!/usr/bin/env python3
"""Propagation over F_p for the whole atom system.

Justification: all 3707 reduction handles are single-use (they occur in exactly
one atom each), so `X - p*H = 0` is exactly the assertion `X == 0 (mod p)` with
H a free quotient.  Every other atom is a polynomial relation that holds over Z
and hence mod p.  So a mod-p solution is a necessary condition, and (modulo the
integer lift) a sufficient one.
"""
import pickle, os, collections, sys, random
from model import Model

HERE = os.path.dirname(os.path.abspath(__file__))
P = 2**256 - 2**32 - 977
NV = 38748


def sqrt_p(a):
    """Return a square root of a mod p (p % 4 == 3) or None."""
    if a == 0:
        return 0
    r = pow(a, (P + 1) // 4, P)
    return r if r * r % P == a else None


class FpEngine:
    def __init__(self, M):
        self.M = M
        self.var2atoms = collections.defaultdict(list)
        for i, vs in enumerate(M.avars):
            for x in vs:
                self.var2atoms[x].append(i)
        self.avarlist = [tuple(vs) for vs in M.avars]
        # reduced compiled monomials
        self.comp = []
        for q in M.polys:
            self.comp.append(tuple((m, c % P) for m, c in q.items() if c % P))
        # handles: vars that occur in exactly one atom, whose coefficient is 0 mod p
        self.handles = set()
        for v, ats in self.var2atoms.items():
            if len(ats) == 1:
                a = ats[0]
                # v's coefficient mod p vanishes?
                tot = [c for m, c in self.comp[a] if v in m]
                if not tot:
                    self.handles.add(v)

    def eval_atom(self, a, val):
        s = 0
        for m, c in self.comp[a]:
            t = c
            for x in m:
                t = t * val[x] % P
            s += t
        return s % P

    def reduce(self, a, val, u):
        c2 = c1 = c0 = 0
        for m, c in self.comp[a]:
            k = 0
            t = c
            for x in m:
                if x == u:
                    k += 1
                else:
                    t = t * val[x] % P
            if k == 0:
                c0 += t
            elif k == 1:
                c1 += t
            else:
                c2 += t
        return c2 % P, c1 % P, c0 % P

    def propagate(self, val, Q=None):
        M = self.M
        unk = [0] * M.na
        for a in range(M.na):
            unk[a] = sum(1 for x in self.avarlist[a] if val[x] is None)
        if Q is None:
            Q = collections.deque(a for a in range(M.na) if unk[a] <= 1)
        conflicts = []
        branch = []
        nassign = 0
        while Q:
            a = Q.popleft()
            miss = [x for x in self.avarlist[a] if val[x] is None]
            if len(miss) > 1:
                continue
            if not miss:
                if self.eval_atom(a, val) != 0:
                    conflicts.append(('eval', a))
                continue
            u = miss[0]
            c2, c1, c0 = self.reduce(a, val, u)
            if c2 == 0:
                if c1 == 0:
                    if c0 != 0:
                        conflicts.append(('const', a))
                    continue
                roots = [(-c0) * pow(c1, -1, P) % P]
            else:
                d = (c1 * c1 - 4 * c2 * c0) % P
                s = sqrt_p(d)
                if s is None:
                    conflicts.append(('nonqr', a))
                    continue
                inv = pow(2 * c2 % P, -1, P)
                roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
            if len(roots) > 1:
                branch.append((u, a, roots))
                continue
            val[u] = roots[0]
            nassign += 1
            for b in self.var2atoms[u]:
                unk[b] -= 1
                if unk[b] <= 1:
                    Q.append(b)
        return nassign, conflicts, branch


def main():
    M = Model()
    E = FpEngine(M)
    print("handles (single-occurrence, p-killed):", len(E.handles))
    val = [None] * NV
    for h in E.handles:
        val[h] = 0     # handles are free mod p; pin them to 0 for now
    rounds = 0
    while True:
        n, conf, br = E.propagate(val)
        rounds += 1
        known = sum(1 for x in val if x is not None)
        print(f"round {rounds}: +{n} known={known} conf={len(conf)} branch={len(br)}", flush=True)
        if conf:
            cc = collections.Counter(k for k, _ in conf)
            print("  conflict kinds:", cc)
            for k, a in conf[:8]:
                print("   ", k, a, M.src[a][:120])
            break
        if not br:
            break
        # branch policy: prefer 0/1 roots (boolean atoms) -> pick 0
        ch = 0
        for u, a, roots in br:
            if val[u] is not None:
                continue
            val[u] = roots[0]
            ch += 1
        print(f"   branch-assigned {ch}")
        if ch == 0:
            break
    known = sum(1 for x in val if x is not None)
    print(f"FINAL mod-p known {known}/{NV}")
    pickle.dump(val, open(os.path.join(HERE, 'fp0.pkl'), 'wb'))


if __name__ == '__main__':
    main()
