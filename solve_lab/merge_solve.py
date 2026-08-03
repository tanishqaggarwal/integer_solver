#!/usr/bin/env python3
"""Solve past the trapdoor by MERGING the twist variables and running strong
propagation mod P.

Merge x_18274->x_9770, x_17728->x_3183 (bakes in the twist), set the 1103 pins,
then propagate to fixpoint with:
  - unit propagation (atom -> 1 unknown, linear or solvable)
  - BILINEAR (product atom with 2 knowns -> the 3rd, via division mod P)
  - residue read-off (bit*(x_B-H)=s*x_C : any 2 known -> the 3rd / the bit)
Report how many vars/bits get determined and whether it's consistent. If it
determines the control bits, verify in Z.
"""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv

MERGE = {18274: 9770, 17728: 3183}   # var -> representative

def remap(v): return MERGE.get(v, v)

def main():
    t0 = time.time()
    atoms0 = load_atoms()
    control = set(json.load(open('control_bits.json')))
    # remap atoms
    atoms = []
    for poly in atoms0:
        out = defaultdict(int)
        for m, c in poly.items():
            nm = tuple(sorted(remap(x) for x in m))
            out[nm] = (out[nm] + c) % P
        atoms.append({m: c % P for m, c in out.items() if c % P})
    print(f"remapped {len(atoms)} atoms ({time.time()-t0:.0f}s)", flush=True)

    val = [None]*NVARS
    var_atoms = defaultdict(list)
    for a, poly in enumerate(atoms):
        for v in set().union(*[set(m) for m in poly]) if poly else []:
            var_atoms[v].append(a)

    contra = [0]
    wl = deque(range(len(atoms))); inwl = [True]*len(atoms)
    def assign(v, x):
        x %= P
        if val[v] is not None:
            if val[v] != x: contra[0] += 1
            return
        val[v] = x
        for a in var_atoms[v]:
            if not inwl[a]: inwl[a] = True; wl.append(a)

    def subst(poly):
        out = defaultdict(int)
        for m, c in poly.items():
            cc = c; nm = []
            for v in m:
                if val[v] is not None: cc = (cc*val[v]) % P
                else: nm.append(v)
            out[tuple(sorted(nm))] = (out[tuple(sorted(nm))] + cc) % P
        return {m: c for m, c in out.items() if c % P}

    # pins: atoms of form (x)-1 or 1*x = c
    steps = 0
    while wl:
        a = wl.popleft(); inwl[a] = False; steps += 1
        r = subst(atoms[a])
        uv = set().union(*[set(m) for m in r]) if r else set()
        if len(uv) == 0:
            if r.get((), 0) % P: contra[0] += 1
            continue
        if len(uv) == 1:
            u = next(iter(uv))
            c0 = c1 = c2 = 0
            for m, c in r.items():
                if len(m) == 0: c0 = (c0+c) % P
                elif len(m) == 1: c1 = (c1+c) % P
                else: c2 = (c2+c) % P
            if c2 == 0 and c1:
                assign(u, (-c0*inv(c1)) % P)
            elif c2 == 0 and c1 == 0:
                if c0 % P: contra[0] += 1
            elif c0 == 0 and (c1+c2) % P == 0:
                pass  # x^2=x boolean -> leave for branching (do NOT guess)
            elif c2 and c1 == 0 and c0:
                # c2*x^2 + c0 = 0 -> x^2 = -c0/c2 ; only if a square root exists (skip; rare)
                pass
            # else quadratic: skip
        elif len(uv) == 2:
            # bilinear: c*x*y + d*x + e*y + f = 0 with one of x,y known handled above;
            # here both unknown -> only solve if it's a pure product a*x*y = b*z form after subst,
            # i.e. monomials are {(x,y)} and {(x,)} etc. Try: if r has a single degree-2 term x*y
            # and the rest are linear in x and y, can't uniquely solve. Skip (needs branching).
            pass
    ndet = sum(1 for x in val if x is not None)
    ndetbits = sum(1 for b in control if val[b] is not None)
    print(f"propagation: {steps} steps, determined {ndet}/{NVARS} vars, {ndetbits}/{len(control)} control bits, contradictions {contra[0]} ({time.time()-t0:.0f}s)", flush=True)
    # how many atoms still unsatisfied (nonzero residual with current partial, zero-filling)?
    v0 = [x if x is not None else 0 for x in val]
    viol = 0
    for poly in atoms:
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t = (t*v0[x]) % P
            s = (s+t) % P
        if s % P: viol += 1
    print(f"with zero-fill: {viol} atoms violated mod P ({time.time()-t0:.0f}s)", flush=True)
    # save determined bits
    detbits = {b: val[b] for b in control if val[b] is not None}
    print(f"determined control bits: {detbits if len(detbits)<40 else len(detbits)}", flush=True)

if __name__ == '__main__':
    main()
