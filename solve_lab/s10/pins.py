"""S10 step 3: for every variable feeding a nonzero atom, show ALL atoms it occurs in,
whether each is currently zero, and whether that atom is 'live' (would actually change)
if we perturbed the variable -- i.e. the effective, not syntactic, constraint set."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
NZ = {22229, 22230, 35758, 35759, 35760, 35761, 35762}

TARGETS = [9118, 8731, 642, 28730, 7068, 2099, 17325, 9413, 1329, 10903, 29854, 31864, 7075]

def live_coeff(a, u):
    """d(atom a)/d(x_u) evaluated at v -- exact, since atoms are multilinear-ish deg<=2."""
    s = 0
    for m, c in L.polys[a].items():
        if u not in m: continue
        t = c
        cnt = 0
        for w in m:
            if w == u and cnt == 0:
                cnt = 1
                continue
            t *= v[w]
        if m.count(u) == 2:      # u^2 term -> derivative 2*c*u
            t = 2 * c * v[u]
        s += t
    return s

for u in TARGETS:
    print(f'\n=== x_{u}  val={v[u]}  ({"FREE" if u not in L.definer else "gate def by atom %d"%L.definer[u]})')
    print(f'    val mod p = {v[u] % P}')
    for a in sorted(L.var_atoms[u]):
        lc = live_coeff(a, u)
        tag = 'NONZERO-ATOM' if a in NZ else ('zero' if av[a] == 0 else 'nz(other)')
        d = 'DORMANT (d/dx=0)' if lc == 0 else f'live d/dx={str(lc)[:40]}'
        own = ' [its definer]' if L.definer.get(u) == a else ''
        print(f'    atom {a:>6} {tag:<13} {d:<52} neqs={len(L.atom2eq.get(a,{}))}{own}')
        print(f'        {L.atom_src[a][:130]}')
