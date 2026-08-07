#!/usr/bin/env python3
"""Exact integer propagation engine over the atom system.

State: val[v] = int or None.  Domain restriction for booleans.
Rules: any atom with <=1 unknown variable is solved exactly (linear or quadratic
over Z).  Conflicts are reported with the atom that caused them.
"""
import pickle, os, collections, math, sys
from model import Model

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


def isqrt_exact(n):
    if n < 0:
        return None
    r = math.isqrt(n)
    return r if r * r == n else None


class Engine:
    def __init__(self, M):
        self.M = M
        self.var2atoms = collections.defaultdict(list)
        for i, vs in enumerate(M.avars):
            for x in vs:
                self.var2atoms[x].append(i)
        # bucket atoms by variable-count
        self.avarlist = [tuple(vs) for vs in M.avars]

    def reduce(self, a, val, u):
        """Return (c2,c1,c0) of atom a in the single unknown u."""
        c2 = c1 = c0 = 0
        for m, c in self.M.compiled[a]:
            k = 0
            t = c
            for x in m:
                if x == u:
                    k += 1
                else:
                    t *= val[x]
            if k == 0:
                c0 += t
            elif k == 1:
                c1 += t
            else:
                c2 += t
        return c2, c1, c0

    def propagate(self, val, queue=None, verbose=False):
        """val: list of int|None.  Returns (status, info)."""
        M = self.M
        unk = [0] * M.na
        for a in range(M.na):
            unk[a] = sum(1 for x in self.avarlist[a] if val[x] is None)
        Q = collections.deque(a for a in range(M.na) if unk[a] <= 1)
        conflicts = []
        nassign = 0
        branchable = []
        while Q:
            a = Q.popleft()
            miss = [x for x in self.avarlist[a] if val[x] is None]
            if len(miss) > 1:
                continue
            if not miss:
                if M.atom_val(a, val) != 0:
                    conflicts.append(('eval', a))
                continue
            u = miss[0]
            c2, c1, c0 = self.reduce(a, val, u)
            roots = None
            if c2 == 0:
                if c1 == 0:
                    if c0 != 0:
                        conflicts.append(('const', a))
                    continue
                if c0 % c1:
                    conflicts.append(('nondiv', a))
                    continue
                roots = [-c0 // c1]
            else:
                d = c1 * c1 - 4 * c2 * c0
                s = isqrt_exact(d)
                if s is None:
                    conflicts.append(('noroot', a))
                    continue
                rs = set()
                for sg in (s, -s):
                    num = -c1 + sg
                    if num % (2 * c2) == 0:
                        rs.add(num // (2 * c2))
                if not rs:
                    conflicts.append(('nonintroot', a))
                    continue
                roots = sorted(rs)
            if len(roots) > 1:
                branchable.append((u, a, roots))
                continue
            val[u] = roots[0]
            nassign += 1
            for b in self.var2atoms[u]:
                unk[b] -= 1
                if unk[b] <= 1:
                    Q.append(b)
        return nassign, conflicts, branchable


def main():
    M = Model()
    E = Engine(M)
    val = [None] * NV
    n, conf, br = E.propagate(val)
    known = sum(1 for x in val if x is not None)
    print(f"assigned {n}, known {known}/{NV}, conflicts {len(conf)}, branchable {len(br)}")
    cc = collections.Counter(k for k, _ in conf)
    print("conflict kinds:", cc)
    for k, a in conf[:15]:
        print("  ", k, a, M.src[a][:120])
    pickle.dump(val, open(os.path.join(HERE, 'prop0.pkl'), 'wb'))


if __name__ == '__main__':
    main()
