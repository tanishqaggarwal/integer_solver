"""CL step 3b: the broadcast class of the 296-bit constants; move the whole class."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
def score(v): return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))
S0 = score(v0)

byval = collections.defaultdict(list)
for u in range(L.NVARS):
    byval[v0[u]].append(u)

for rhs in (14623, 14853):
    val = v0[rhs]
    grp = byval[val]
    fg = [u for u in grp if u in FREE]
    print(f'\nx_{rhs}: value has {len(grp)} variables sharing it, {len(fg)} of them FREE')
    print(f'   free members: {fg}')
    print(f'   computed members: {[u for u in grp if u not in FREE][:30]}')

# how many distinct 296-bit values are there?
big = {k: g for k, g in byval.items() if k and abs(k).bit_length() > 200 and abs(k).bit_length()<400}
print(f'\ndistinct values in 200..400 bits: {len(big)}; class sizes: '
      f'{sorted((len(g) for g in big.values()), reverse=True)[:20]}')
