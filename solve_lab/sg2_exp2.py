#!/usr/bin/env python3
"""Experiment 2: after the residue-fix, find which ATOMS are nonzero (the genuine ripple checks)."""
import heal_harness as H
import sg2_lib as L
p = H.p
atoms = L.load_atoms_full()
A = {a['idx']: a for a in atoms}

def evatom(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for vv in m: t *= val[vv]
        s += t
    return s

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()

x2099 = H.val[2099]; x19964 = H.val[19964]
H.val[7068] = x2099
H.val[4432] = x19964
H.val[17325] = 0; H.val[9413] = 0
H.forward()

nz = []
for a in atoms:
    r = evatom(a['poly'], H.val)
    if r != 0:
        nz.append((a['idx'], r))
print(f"nonzero atoms after residue-fix: {len(nz)}")
for ai, r in nz:
    a = A[ai]
    deg = max((len(m) for m in a['poly']), default=0)
    print(f"  atom {ai} (deg{deg}, n_eq={a['n_eq']}, eqs={a['eqs']}): {a['repr'][:110]}")
    print(f"      value/p={r//p}  value%p={r%p}")
