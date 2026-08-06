import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, BOOLATOM, CANON, F2, pot, FORBID
P = 2**256-2**32-977

w26 = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w26)
v0  = L.load(os.path.join(HERE,'mod9118_0.json')); CANON.fwd(v0)

SEVEN=[22229,22230,35758,35759,35760,35761,35762]
print('=== FRAME2 cone of the seven ===')
c2 = F2.cone(SEVEN)
print(f'{len(c2)} vars: {sorted(c2)}')
print(f'free in F2: {sorted(c2 & F2.FREE)}')
print(f'booleans: {sorted(c2 & BOOL)}')
print(f'FORBID inside: {sorted(c2 & FORBID)}')
print('values:')
for u in sorted(c2):
    d = F2.definer.get(u)
    src = L.atom_src[d][:90] if d is not None else 'FREE'
    print(f'  x_{u:<6} = {str(w26[u])[:34]:<34} bits={w26[u].bit_length():<4} {src}')

print('\n=== CANONICAL cone of a21617,a29539 : boolean free inputs ===')
cc = CANON.cone([21617,29539])
bc = sorted(cc & BOOL & CANON.FREE - FORBID)
print(f'{len(cc)} vars, {len(bc)} boolean free inputs')
print(f'their current values in mod9118_0: {[(u, v0[u]) for u in bc]}')
ones = [u for u in bc if v0[u]==1]; zers=[u for u in bc if v0[u]==0]
print(f'  ones {len(ones)}: {ones}')
print(f'  zeros {len(zers)}: {zers}')
