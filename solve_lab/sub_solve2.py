#!/usr/bin/env python3
"""Enumerate the 22-bit subproblem controlling all 4 violated atoms.
Fast mod-P sub-circuit eval (precomputed inverses); collect patterns zeroing all
4 targets; verify each hit over Z with full propagation + exact checker."""
import json, time, sys
from collections import deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import boolean_vars, ProvEngine
from modp import ModPEngine, P, NVARS, inv

BITS22 = [1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,
          27512,29682,30104,30596,30658,30792,33251,37748,37885,38116]
TARGETS = [27973, 27978, 41470, 45004]

def main():
    t0 = time.time()
    atoms = load_atoms(); bset = boolean_vars(atoms)
    base = ModPEngine(atoms); base.propagate()
    for b in [x for x in bset if base.val[x] is None]:
        if base.val[b] is None: base.assign(b, 0); base.propagate()
    for v in range(NVARS):
        if base.val[v] is None: base.assign(v, 0); base.propagate()
    base_val = [x if x is not None else 0 for x in base.val]

    e2 = ModPEngine(atoms)
    for b in BITS22: e2.assign(b, 1)
    e2.propagate()
    for b in [x for x in bset if e2.val[b] is None]:
        if e2.val[b] is None: e2.assign(b, 0); e2.propagate()
    for v in range(NVARS):
        if e2.val[v] is None: e2.assign(v, 0); e2.propagate()
    allon = [x if x is not None else 0 for x in e2.val]
    cone = set(v for v in range(NVARS) if allon[v] != base_val[v])

    eo = json.load(open('eval_order.json')); order = eo['order']; prov = eo['prov']
    coneb = cone | set(BITS22)
    suborder = []
    defatom = {}
    for v in order:
        if v in coneb and prov[v] and prov[v][0] >= 0:
            suborder.append(v)
            defatom[v] = [(c % P, tuple(m)) for m, c in atoms[prov[v][0]].items()]
    tgt = {ai: [(c % P, tuple(m)) for m, c in atoms[ai].items()] for ai in TARGETS}
    # precompute constant inverse for wires whose target-term is just (v,)
    neg_inv = {}
    for v in suborder:
        tt = [(coef, m) for coef, m in defatom[v] if v in m]
        if len(tt) == 1 and tt[0][1] == (v,):
            neg_inv[v] = (-inv(tt[0][0])) % P
    print(f"cone {len(cone)}, suborder {len(suborder)}, const-inv wires {len(neg_inv)} ({time.time()-t0:.0f}s)", flush=True)

    val = list(base_val)
    def evalpat_targets(code):
        for i, b in enumerate(BITS22): val[b] = (code >> i) & 1
        for v in suborder:
            da = defatom[v]; ni = neg_inv.get(v)
            if ni is not None:
                c0 = 0
                for coef, m in da:
                    if v in m: continue
                    t = coef
                    for x in m: t = (t*val[x]) % P
                    c0 = (c0 + t) % P
                val[v] = (c0*ni) % P
            else:
                c0 = c1 = 0
                for coef, m in da:
                    if v in m:
                        t = coef
                        for x in m:
                            if x != v: t = (t*val[x]) % P
                        c1 = (c1 + t) % P
                    else:
                        t = coef
                        for x in m: t = (t*val[x]) % P
                        c0 = (c0 + t) % P
                val[v] = (-c0*inv(c1)) % P if c1 else 0
        for ai in TARGETS:
            s = 0
            for coef, m in tgt[ai]:
                t = coef
                for x in m: t = (t*val[x]) % P
                s = (s+t) % P
            if s: return False
        return True

    # sanity check vs propagation
    print(f"sanity all-on: {'all4=0' if evalpat_targets((1<<22)-1) else 'not all 0 (matches prop: 45004 breaks)'}", flush=True)
    t1 = time.time(); hits = []; N = 1 << 22
    for code in range(N):
        if evalpat_targets(code):
            hits.append(code)
            ones = [b for i, b in enumerate(BITS22) if (code >> i) & 1]
            print(f"  HIT {code:#x}: bits {ones}", flush=True)
        if code and code % 500000 == 0:
            print(f"  ...{code}/{N} ({time.time()-t1:.0f}s) hits={len(hits)}", flush=True)
    print(f"enum done {time.time()-t1:.0f}s; {len(hits)} hits", flush=True)
    json.dump([[b for i,b in enumerate(BITS22) if (h>>i)&1] for h in hits], open('sub2_hits.json','w'))

    # verify hits over Z
    A = load_atoms()
    for h in hits[:200]:
        ones = [b for i, b in enumerate(BITS22) if (h >> i) & 1]
        eng = ProvEngine(A); eng.propagate()
        for b in ones: eng.assign(b, 1, ('f', ()))
        eng.propagate()
        for b in [x for x in bset if eng.val[x] is None]:
            if eng.val[b] is None: eng.assign(b, 0, ('f', ())); eng.propagate()
        for v in range(NVARS):
            if eng.val[v] is None: eng.assign(v, 0, ('f', ())); eng.propagate()
        v = [x if x is not None else 0 for x in eng.val]
        vio = 0
        for poly in A:
            s = 0
            for m, c in poly.items():
                t = c
                for x in m: t *= v[x]
                s += t
            if s: vio += 1
        print(f"  Z-verify bits {ones}: {vio} atoms violated", flush=True)
        if vio == 0:
            json.dump({f"x_{i}": v[i] for i in range(NVARS)}, open('cand_SUBSOLVE.json', 'w'))
            print("  *** FULL SOLUTION over Z! wrote cand_SUBSOLVE.json ***", flush=True)
            return

if __name__ == '__main__':
    main()
