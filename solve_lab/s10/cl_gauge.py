"""CL: global gauge shift of the broadcast constant classes K1 / K2."""
import os, sys, json, collections, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)
print(f'base score {S0}  nz {nz0}')

K1FREE = [8778, 14623, 16742, 31339, 33462]
K2FREE = [14853, 22152, 22649]
tgt1 = v0[24548] % P            # what a21617/a7930 want K1 to be
d1 = (tgt1 - v0[14623]) % P
print(f'K1 = {v0[14623]%P}\nwant {tgt1}\ndelta1 = {d1}')

def run(shift, tag, absorb_all=True):
    v = list(v0)
    for u, d in shift.items(): v[u] = v[u] + d
    ad.fwd(v, rounds=10)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    # absorb every nonzero CHECK with a p-handle
    if absorb_all:
        for _ in range(3):
            av = L.all_atom_values(v)
            nz = [a for a in range(L.NA) if av[a] and a not in atom_out]
            prog = False
            for a in nz:
                if av[a] % P == 0 or True:
                    if E.absorb(v, a): prog = True
            if not prog: break
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    bad = set(L.failing_eqs(av))
    s = L.NEQ - len(bad)
    print(f'{tag}: score {s} ({s-S0:+d})  nz({len(nz)})={nz[:14]}  broke {len(bad-bad0)} fixed {len(bad0-bad)}')
    return v, s, nz

# 1. shift only x_14623
run({14623: d1}, 'x_14623 alone')
# 2. shift all K1 free members
run({u: d1 for u in K1FREE}, 'ALL K1 free members')
# 3. subsets
for k in range(1, 5):
    for sub in itertools.combinations(K1FREE, k):
        if 14623 not in sub: continue
        run({u: d1 for u in sub}, f'K1 subset {sub}')
