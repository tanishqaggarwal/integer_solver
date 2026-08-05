#!/usr/bin/env python3
"""Iterative forward-model healer starting from agentA.
Fix G1,G2 residues; then repeatedly heal nonzero deg-1 check atoms using free-input slacks.
Log the propagation chain; detect termination / cycle / wall."""
import heal_harness as H
import sg2_lib as L
import pickle, sys
p = H.p
atoms = L.load_atoms_full()
A = {a['idx']: a for a in atoms}
idx = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/idx.pkl','rb'))
var_atoms = idx['var_atoms']
natoms = {v: len(var_atoms[v]) for v in var_atoms}

def evatom(poly, val):
    s = 0
    for m, c in poly.items():
        t = c
        for vv in m: t *= val[vv]
        s += t
    return s

def deg(poly):
    return max((len(m) for m in poly), default=0)

def nonzero_atoms():
    out = []
    for a in atoms:
        r = evatom(a['poly'], H.val)
        if r != 0:
            out.append((a['idx'], r))
    return out

# --- linear coefficient of a single free var v in atom (after forward, treating others fixed) ---
def lin_coef(poly, v):
    """coefficient of v treating current val; returns (c1, c2) for c2*v^2+c1*v+... ; but v is free."""
    c1 = 0; c2 = 0
    for m, c in poly.items():
        cnt = m.count(v)
        if cnt == 0: continue
        rest = [u for u in m if u != v]
        prod = c
        for u in rest: prod *= H.val[u]
        if cnt == 1: c1 += prod
        elif cnt == 2: c2 += prod
    return c1, c2

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()
# fix residues
H.val[7068] = H.val[2099]; H.val[4432] = H.val[19964]
H.val[17325] = 0; H.val[9413] = 0
H.forward()

used = set()  # free vars we've moved
log = []
for it in range(200):
    nz = nonzero_atoms()
    if not nz:
        print(f"iter {it}: ALL ATOMS ZERO! SOLVED core.")
        break
    deg1 = [(ai, r) for ai, r in nz if deg(A[ai]['poly']) == 1]
    deg2 = [(ai, r) for ai, r in nz if deg(A[ai]['poly']) == 2]
    # pick a deg-1 check to heal: find its best free slack (unit coef, fewest atoms, not yet a dead-end)
    healed = False
    # sort deg1 by ... just take first
    for ai, r in deg1:
        poly = A[ai]['poly']
        frees = [v for v in L.atom_vars(poly) if v in H.freeinp]
        # candidate slacks: free var with |c1|==1, c2==0
        cands = []
        for v in frees:
            c1, c2 = lin_coef(poly, v)
            if c2 == 0 and abs(c1) == 1:
                cands.append((natoms.get(v, 99), v, c1))
        cands.sort()
        for _, v, c1 in cands:
            # set v to zero the atom: current r = c1*val[v] + rest ; rest = r - c1*val[v]
            # want c1*newv + rest = 0 -> newv = -rest/c1 = val[v] - r/c1
            if r % c1 != 0: continue
            newv = H.val[v] - r // c1
            H.val[v] = newv
            H.forward()
            log.append((it, ai, v, natoms.get(v)))
            healed = True
            break
        if healed: break
    if not healed:
        print(f"iter {it}: {len(nz)} nonzero, deg1={len(deg1)} deg2={len(deg2)}; NO healable deg-1 slack.")
        print(f"   deg1 atoms: {[ai for ai,_ in deg1][:10]}")
        print(f"   deg2 atoms: {[ai for ai,_ in deg2][:10]}")
        break
    if it % 1 == 0:
        print(f"iter {it}: healed atom {log[-1][1]} via free x_{log[-1][2]} ({log[-1][3]} atoms); now {len(nonzero_atoms())} nonzero")

print(f"\nchain length: {len(log)}")
print("free vars used (in order):", [v for _,_,v,_ in log][:60])
