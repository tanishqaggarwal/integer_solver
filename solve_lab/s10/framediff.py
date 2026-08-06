import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from newton import score
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
A = L.load(os.path.join(HERE, 'mod9118_0.json'))          # canonical, 39009
B = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
df = [u for u in range(L.NVARS) if A[u] != B[u]]
dfree = [u for u in df if u in FREE]
print(f'variables differing: {len(df)};  of which FREE INPUTS: {len(dfree)}')
print(f'free-input differences: {dfree}')
for u in dfree:
    print(f'  x_{u:<7} canon={str(A[u])[:26]:<28} delivered={str(B[u])[:26]:<28} '
          f'consumers {len(L.var_atoms[u])}  equal mod p? {A[u] % P == B[u] % P}')
# what happens if we import the delivered free inputs into the canonical frame?
v = list(A)
for u in dfree: v[u] = B[u]
ad.fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'\nimport ALL delivered free inputs -> canonical: score {score(v)} '
      f'nonzero {[a for a in range(L.NA) if av[a]]}')
# one at a time
print('\none free input at a time:')
for u in dfree:
    w = list(A); w[u] = B[u]
    ad.fwd(w, rounds=8)
    aw = L.all_atom_values(w)
    print(f'  x_{u:<7} -> score {score(w):>6}  nonzero {[a for a in range(L.NA) if aw[a]]}')
