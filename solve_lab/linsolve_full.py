#!/usr/bin/env python3
"""Full sparse Gaussian elimination over GF(p) on ALL linear atoms (degree-1).
Determines the TRUE linear rank and the free-variable set (stronger than unit
propagation, which only handles 1-unknown atoms). Reports whether the twist vars
x_9770,x_18274,x_3183,x_17728 are linearly determined by the free set, and how
many free vars remain -> the true dimension of the nonlinear core."""
import json, time
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS

P = (1 << 61) - 1
def inv(a): return pow(a % P, P-2, P)

def main():
    t0 = time.time()
    atoms = load_atoms()
    # collect purely-linear atoms (all monomials degree<=1)
    lin = []
    for poly in atoms:
        if all(len(m) <= 1 for m in poly):
            row = {}
            const = 0
            for m, c in poly.items():
                if len(m) == 0: const = (const + c) % P
                else: row[m[0]] = (row.get(m[0], 0) + c) % P
            row = {v: c for v, c in row.items() if c % P}
            if row: lin.append((row, (-const) % P))
    print(f"linear atoms: {len(lin)} ({time.time()-t0:.0f}s)", flush=True)

    # sparse Gaussian elimination over GF(P)
    pivot_of = {}       # var -> (row dict, rhs) with that var as pivot (coef 1)
    rank = 0
    incons = 0
    for row, rhs in lin:
        r = dict(row); rr = rhs
        # reduce by known pivots
        changed = True
        while changed:
            changed = False
            for v in list(r):
                if v in pivot_of:
                    f = r.pop(v); prow, prhs = pivot_of[v]
                    for k, cv in prow.items():
                        nv = (r.get(k, 0) - f*cv) % P
                        if nv: r[k] = nv
                        elif k in r: del r[k]
                    rr = (rr - f*prhs) % P
                    changed = True
                    break
        if not r:
            if rr % P: incons += 1
            continue
        # pick pivot = var with smallest index
        pv = min(r); f = inv(r[pv])
        prow = {k: (cv*f) % P for k, cv in r.items() if k != pv}
        pivot_of[pv] = (prow, (rr*f) % P)
        rank += 1
        if rank % 2000 == 0: print(f"  rank {rank} ({time.time()-t0:.0f}s)", flush=True)
    print(f"linear RANK {rank}, inconsistent {incons} ({time.time()-t0:.0f}s)", flush=True)

    # variables appearing in linear atoms
    linvars = set()
    for row, rhs in lin: linvars.update(row)
    pivoted = set(pivot_of)
    freelin = linvars - pivoted
    print(f"vars in linear atoms: {len(linvars)}; pivoted (determined): {len(pivoted)}; free-in-linear: {len(freelin)}", flush=True)
    print(f"total vars NOT in any linear atom: {NVARS - len(linvars)} (free too)", flush=True)
    print(f"=> TRUE free dimension (linear) approx: {NVARS - rank}", flush=True)

    for w in (9770, 3183, 18274, 17728, 8821, 6773, 17233):
        if w in pivot_of:
            prow, prhs = pivot_of[w]
            nz = len(prow)
            controlvars = set(json.load(open('control_bits.json')))
            deps_ctrl = [v for v in prow if v in controlvars]
            print(f"  x_{w}: linearly DETERMINED by {nz} vars (rhs={prhs}); control-bit deps: {len(deps_ctrl)}")
        elif w in linvars:
            print(f"  x_{w}: FREE (appears in linear atoms but is a pivot-free var)")
        else:
            print(f"  x_{w}: not in any linear atom (free / nonlinear-only)")
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
