#!/usr/bin/env python3
"""Greedy heal from agentA: fix G1,G2 residues, then heal ripple checks using free slacks."""
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

def nonzero_atoms():
    return [(a['idx'], evatom(a['poly'], H.val)) for a in atoms if evatom(a['poly'], H.val) != 0]

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()

# fix G1,G2
H.val[7068] = H.val[2099]; H.val[4432] = H.val[19964]
H.val[17325] = 0; H.val[9413] = 0
H.forward()
nz = nonzero_atoms()
print(f"after residue-fix: {len(nz)} nonzero atoms: {[ai for ai,_ in nz]}")

# Try setting x_2964 (free) to fix atom 7450, x_24548 (free) to fix atom 7452
# atom 7450: x_2964 - x_26756 - x_579  -> x_2964 = x_26756 + x_579
H.val[2964] = H.val[26756] + H.val[579]
# atom 7452: 9367949*(x_24548 - x_25442) - x_7927 -> x_24548 = x_25442 + x_7927/9367949
num = H.val[7927]
if num % 9367949 == 0:
    H.val[24548] = H.val[25442] + num // 9367949
    print(f"set x_24548 exactly")
else:
    print(f"x_7927={num} not divisible by 9367949 (rem {num%9367949}); need x_7927 adjustable")
H.forward()
nz = nonzero_atoms()
print(f"after setting x_2964,x_24548: {len(nz)} nonzero: {[(ai, r//p, r%p) for ai,r in nz][:20]}")

# report x_2964 and x_24548 roles
import pickle
idx = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/idx.pkl','rb'))
var_atoms = idx['var_atoms']
for v in [2964, 24548, 7927, 11052, 15616]:
    kind='FREE' if v in H.freeinp else 'gate'
    print(f"x_{v} [{kind}] in atoms: {var_atoms[v]}")
