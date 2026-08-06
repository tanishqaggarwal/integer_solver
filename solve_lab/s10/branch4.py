"""S10 step 91: the (x_2081, x_4287) MUX has FOUR branches. We have only ever
been in (1,0).  Branch (1,1) sets x_21279 = 1, hence x_7075 = 0, which kills the
multiplier in a35759 and a35761 -- and through them a35758 and a35760.

Test each branch from the CANONICAL frame with full forward re-evaluation.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
base = L.load(os.path.join(HERE, 'forward_state.json'))
av0 = L.all_atom_values(base)
f0 = L.failing_eqs(av0)
nz0 = [a for a in range(L.NA) if av0[a]]
print(f'canonical frame: failing {len(f0)}, nonzero atoms {len(nz0)}')
print(f'  nonzero: {nz0}')
print(f'  controls now: x_2081={base[2081]} x_4287={base[4287]} '
      f'x_21279={base[21279]} x_7075={base[7075]}')

for b1 in (0, 1):
    for b2 in (0, 1):
        v = list(base)
        v[2081] = b1; v[4287] = b2
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        fail = L.failing_eqs(av)
        nz = [a for a in range(L.NA) if av[a]]
        print(f'\n(x_2081,x_4287)=({b1},{b2}): x_21279={v[21279]} x_7075={v[7075]} '
              f'x_2099={str(v[2099])[:20]}')
        print(f'   failing eqs {len(fail):>5}   nonzero atoms {len(nz):>4}   '
              f'score {L.NEQ - len(fail)}')
        print(f'   nonzero atoms: {nz[:40]}')
