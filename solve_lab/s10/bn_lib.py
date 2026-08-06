"""bn_lib: shared helpers for the boolean-carrier sweep.

The 39026 partial is OFF-MANIFOLD: 5 gate atoms are deliberately nonzero.
Plain ad.fwd repairs them and loses 30 equations.  fwdb() forward-evaluates
while leaving the blocked gates alone.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

BEST = os.path.join(LAB,'best','new_instance_partial_39026.json')
P = 2**256 - 2**32 - 977
FREESET = set(ad.FREE)

def broken_gates(v=None):
    if v is None: v = L.load(BEST)
    av = L.all_atom_values(v)
    return [a for a in range(L.NA) if av[a] and a in L.atom_out]

def fwdb(v, block, rounds=1):
    """forward evaluate skipping variables defined by a blocked gate atom"""
    blocked_vars = set()
    for a in block:
        if a in L.atom_out: blocked_vars.add(L.atom_out[a][1])
    for _ in range(rounds):
        for u in ad.ORDER:
            if u in blocked_vars: continue
            nv = T.solve_lin(L.definer[u], u, v)
            if nv is not None: v[u] = nv
    return v

def score(v):
    av = L.all_atom_values(v)
    f = L.failing_eqs(av)
    return L.NEQ-len(f), f, av

def bools_map():
    cen = json.load(open(os.path.join(HERE,'bn_census.json')))
    return {int(a):tuple(t) for a,t in cen['bools'].items()}
