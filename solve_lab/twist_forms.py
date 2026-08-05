#!/usr/bin/env python3
"""Fully reduce x_9770,x_3183,x_18274,x_17728 to free vars via the linear system,
then trace each contributing var: is it a residue-load output (bit*(x_B-HUGE)=s*x_C)?
If the twist reduces to a relation among residue loads, it becomes a subset-sum in
the bits (solvable). Also express the twist gaps in the reduced basis."""
import json, time
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
from confluent_eval5 import boolean_vars

P = (1 << 61) - 1
def inv(a): return pow(a % P, P-2, P)

def build_pivots(atoms):
    lin = []
    for poly in atoms:
        if all(len(m) <= 1 for m in poly):
            row = {}; const = 0
            for m, c in poly.items():
                if len(m) == 0: const = (const+c) % P
                else: row[m[0]] = (row.get(m[0], 0)+c) % P
            row = {v: c for v, c in row.items() if c % P}
            if row: lin.append((row, (-const) % P))
    pivot_of = {}
    for row, rhs in lin:
        r = dict(row); rr = rhs; changed = True
        while changed:
            changed = False
            for v in list(r):
                if v in pivot_of:
                    f = r.pop(v); prow, prhs = pivot_of[v]
                    for k, cv in prow.items():
                        nv = (r.get(k, 0) - f*cv) % P
                        if nv: r[k] = nv
                        elif k in r: del r[k]
                    rr = (rr - f*prhs) % P; changed = True; break
        if not r: continue
        pv = min(r); f = inv(r[pv])
        pivot_of[pv] = ({k: (cv*f) % P for k, cv in r.items() if k != pv}, (rr*f) % P)
    return pivot_of

def reduce_var(v, pivot_of, memo, depth=0):
    # returns dict free_var->coef for x_v (fully reduced), plus constant
    if v in memo: return memo[v]
    if v not in pivot_of:
        memo[v] = ({v: 1}, 0); return memo[v]
    prow, prhs = pivot_of[v]
    acc = defaultdict(int); const = prhs
    for k, c in prow.items():
        sub, sc = reduce_var(k, pivot_of, memo, depth+1)
        for fv, fc in sub.items(): acc[fv] = (acc[fv] + c*fc) % P
        const = (const + c*sc) % P
    acc = {k: c for k, c in acc.items() if c % P}
    memo[v] = (acc, const % P); return memo[v]

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    control = set(json.load(open('control_bits.json')))
    pivot_of = build_pivots(atoms)
    print(f"pivots {len(pivot_of)} ({time.time()-t0:.0f}s)", flush=True)
    import sys; sys.setrecursionlimit(100000)
    memo = {}
    forms = {}
    for w in (9770, 3183, 18274, 17728):
        f, c = reduce_var(w, pivot_of, memo)
        forms[w] = (f, c)
        print(f"x_{w} = {len(f)} free vars + const {c}; vars: {sorted(f)[:10]}", flush=True)

    # residue-load map: var -> (bit, HUGE) if it's a simple load output
    load_of = {}
    for a, poly in enumerate(atoms):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2]
        huges = [c for m, c in poly.items() if abs(c) > 10**40]
        if len(prod2) != 1 or not huges: continue
        (m2, c2) = prod2[0]
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in control:
                load_of[xB] = (bit, max(huges, key=abs)); break
    print(f"\nresidue-load outputs mapped: {len(load_of)}", flush=True)

    # for the twist gap x_9770 - x_18274, express in free basis and check if free vars are loads
    def combine(fa, ca, fb, cb):
        acc = defaultdict(int)
        for k, c in fa.items(): acc[k] = (acc[k]+c) % P
        for k, c in fb.items(): acc[k] = (acc[k]-c) % P
        return {k: c for k, c in acc.items() if c % P}, (ca-cb) % P
    for (u, v, name) in [(9770, 18274, 'gap1'), (3183, 17728, 'gap2')]:
        g, gc = combine(forms[u][0], forms[u][1], forms[v][0], forms[v][1])
        loads = [k for k in g if k in load_of]
        nonloads = [k for k in g if k not in load_of]
        print(f"\n{name} = x_{u}-x_{v}: {len(g)} free vars, const {gc}")
        print(f"  free vars that ARE residue loads: {len(loads)}/{len(g)}; non-load free vars: {nonloads[:8]}")
        # which bits gate the load free vars
        bits = sorted(set(load_of[k][0] for k in loads))
        BITS22 = {1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116}
        n22 = sum(1 for b in bits if b in BITS22)
        print(f"  loads gated by {len(bits)} bits ({n22} of the 22-side, {len(bits)-n22} of 233-side)")
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
