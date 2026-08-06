"""S10 step 20: the MUX branch's OWN handles.

Session 9 concluded the MUX branch is closed because the load obligation
(x_4306 == 0, x_27177 == 0 mod p) forces x_9118 to a unique residue that clashes
with the residual congruence.  That analysis treated x_31861 and x_14865 as
constants pinned by the newly-lit load pins a3568/a3570.  But those pins have
handles:  a3568 = x_4287*(x_31861 - C1) - 13479571*x_27676
          a3570 = x_4287*(x_14865 - C2) - x_7574
If x_27676 / x_7574 are free, x_31861 and x_14865 are knobs too -- giving 4 knobs
for 4 conditions instead of 2 for 4.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(HERE, 'mux_on.json'))
av = L.all_atom_values(v)

print('=== handles of the two load pins lit by x_4287 ===')
for a in (3568, 3570):
    print(f'\na{a}: {L.atom_src[a]}')
    for u in sorted(L.avars[a]):
        free = u not in L.definer
        d = L.definer.get(u)
        print(f'   x_{u:<7} free={str(free):<6} natoms={len(L.var_atoms[u]):<3} '
              f'neqs={len(L.var_eqs[u]):<4} val={str(v[u])[:24]}')
        if not free:
            print(f'        definer a{d}: {L.atom_src[d][:110]}')

print('\n=== the load-obligation variables ===')
for u in (4306, 27177, 25490, 37944, 21671, 37530, 9106, 2239, 31731,
          9629, 23754, 35619, 27676, 7574):
    d = L.definer.get(u)
    print(f'  x_{u:<7} free={u not in L.definer} natoms={len(L.var_atoms[u]):<3} '
          f'neqs={len(L.var_eqs[u]):<4} val={str(v[u])[:26]}')
    if d is not None:
        print(f'      a{d}: {L.atom_src[d][:120]}')

print('\n=== sensitivity of the three load atoms to the four candidate knobs ===')
KN = [9118, 8731, 27676, 7574, 31861, 14865]
LOAD = [19088, 22233, 22235]
base = {a: av[a] for a in LOAD}
for u in KN:
    w = list(v)
    L.ripple(w, {u: v[u] + 1})
    wav = L.all_atom_values(w)
    d = {a: (wav[a] - base[a]) for a in LOAD}
    nz = [a for a in range(L.NA) if wav[a]]
    print(f'  x_{u:<7} d(load atoms)={ {a: str(x)[:22] for a,x in d.items()} } '
          f'#nz_atoms={len(nz)}')
