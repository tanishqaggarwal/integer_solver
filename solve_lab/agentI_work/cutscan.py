#!/usr/bin/env python3
"""Cut screen: disable a candidate SET of atoms (allow them nonzero), re-run the
mod-p propagation, and see whether the rest becomes consistent.

A conflict-free run means the whole mod-p defect can be absorbed by the disabled
set; the equation cost is then the number of equations whose atom combination is
nonzero.  Anything with cost < 7 beats the 39,026 deliverable.
"""
import pickle, os, collections, sys, time, json, itertools
from model import Model, load_assign
from fp import FpEngine, P, sqrt_p

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


class Cutter:
    def __init__(self):
        self.M = Model(); self.E = FpEngine(self.M)
        self.wit = load_assign(os.path.join(HERE, '..', 'best',
                                            'new_instance_partial_39026.json'))
        self.witp = [w % P for w in self.wit]
        self.comp = self.E.comp; self.v2a = self.E.var2atoms

    def run(self, disable=(), preassign=None, forcebool=None):
        M = self.M; comp = self.comp; v2a = self.v2a
        dis = bytearray(M.na)
        for a in disable:
            dis[a] = 1
        fb = dict(forcebool or {})
        val = [None] * NV
        reason = [None] * NV
        for k, x in (preassign or {}).items():
            val[k] = x % P; reason[k] = 'pre'
        Q = collections.deque(a for a in range(M.na) if not dis[a])
        inq = bytearray(M.na)
        for a in Q:
            inq[a] = 1
        conflicts = []
        dead = bytearray(M.na)

        def reduce_atom(a):
            red = {}
            for m, c in comp[a]:
                t = c; um = []
                for x in m:
                    if val[x] is None:
                        um.append(x)
                    else:
                        t = t * val[x] % P
                if t == 0:
                    continue
                k = tuple(sorted(um))
                red[k] = (red.get(k, 0) + t) % P
            return {k: c for k, c in red.items() if c}

        def assign(u, x, r):
            val[u] = x; reason[u] = r
            for b in v2a[u]:
                if not inq[b] and not dis[b]:
                    inq[b] = 1; Q.append(b)

        while True:
            branch = []
            while Q:
                a = Q.popleft(); inq[a] = 0
                if dis[a]:
                    continue
                red = reduce_atom(a)
                vs = set()
                for k in red:
                    vs |= set(k)
                if not vs:
                    if red and not dead[a]:
                        dead[a] = 1; conflicts.append(a)
                    continue
                if len(vs) > 1:
                    continue
                u = next(iter(vs))
                c0 = red.get((), 0); c1 = red.get((u,), 0); c2 = red.get((u, u), 0)
                if c2 == 0:
                    assign(u, (-c0) * pow(c1, -1, P) % P, a)
                else:
                    d = (c1 * c1 - 4 * c2 * c0) % P
                    s = sqrt_p(d)
                    if s is None:
                        if not dead[a]:
                            dead[a] = 1; conflicts.append(a)
                        assign(u, (-c1) * pow(2 * c2 % P, -1, P) % P, a)
                        continue
                    inv = pow(2 * c2 % P, -1, P)
                    roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
                    if len(roots) == 1:
                        assign(u, roots[0], a)
                    else:
                        branch.append((u, a, roots))
            pend = [(u, a, r) for u, a, r in branch if val[u] is None]
            if not pend:
                break
            for u, a, roots in pend:
                if val[u] is not None:
                    continue
                x = fb.get(u, self.witp[u])
                assign(u, x if x in roots else roots[0], 'dec')
        return val, conflicts


def main():
    C = Cutter(); M = C.M
    base_val, base_conf = C.run()
    print("baseline conflicts:", base_conf, flush=True)
    # candidate cut atoms: everything in the derivation cone of the conflicts,
    # restricted to atoms appearing in few equations
    mode = sys.argv[1] if len(sys.argv) > 1 else 'single'
    lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hi = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
    maxeq = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    cands = [a for a in range(M.na) if len(M.atom_eqs[a]) <= maxeq]
    print("candidates:", len(cands), flush=True)
    res = {}
    t0 = time.time()
    out = os.path.join(HERE, f'cut_{mode}_{lo}_{hi}_{maxeq}.json')
    for i, a in enumerate(cands[lo:hi]):
        val, conf = C.run(disable=(a,))
        if len(conf) < len(base_conf):
            res[a] = {'nconf': len(conf), 'conf': conf[:10],
                      'eqs': len(M.atom_eqs[a]), 'src': M.src[a][:100]}
            print(f"  HIT a{a} nconf={len(conf)} eqs={len(M.atom_eqs[a])} {M.src[a][:80]}",
                  flush=True)
        if i % 100 == 0:
            json.dump({str(k): v for k, v in res.items()}, open(out, 'w'))
            print(f"  {lo+i}/{min(hi,len(cands))} t={time.time()-t0:.0f}s", flush=True)
    json.dump({str(k): v for k, v in res.items()}, open(out, 'w'))
    print("done", len(res), time.time() - t0)


if __name__ == '__main__':
    main()
