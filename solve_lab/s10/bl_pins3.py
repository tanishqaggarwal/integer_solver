"""bl_pins3: conditional constant LOADS.

Real shape in this instance:   b*(x - K) - c*z = 0
i.e. monomials {(b,x):c1, (b,):-c1*K, (z,):-c}.  With b=1 the wire x is pinned to
K (up to the small slack c*z);  with b=0 the atom degenerates to c*z = 0 and x is
completely UNPINNED.  Enumerate every atom carrying a monomial (b,) whose
coefficient is p-scale, with b boolean.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import BOOL, BOOLATOM, CANON, F2, pot, FORBID
P = 2**256-2**32-977

w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
v0 = L.load(os.path.join(HERE,'mod9118_0.json')); CANON.fwd(v0)

# --- all conditional constant loads ---
LOADS = []          # (atom, b, K, partners, nmon)
for a, p in enumerate(L.polys):
    for m, c in p.items():
        if len(m) == 1 and abs(c) > 2**64:
            b = m[0]
            if b not in BOOL: continue
            partners = sorted(set(z for mm in p for z in mm if b in mm and z != b))
            LOADS.append((a, b, -c, partners, len(p)))
print(f'conditional constant loads (boolean b carrying a p-scale coefficient): {len(LOADS)}')
print(f'  distinct gate booleans: {len(set(t[1] for t in LOADS))}')
print(f'  distinct atoms:         {len(set(t[0] for t in LOADS))}')
byb = collections.Counter(t[1] for t in LOADS)
print(f'  gates: {sorted(byb.items())[:40]}')

clean = [t for t in LOADS if t[4] == 3 and len(t[3]) == 1]
print(f'\nCLEAN pins  b*(x-K) - c*z  (3 monomials, one partner): {len(clean)}')
for a, b, K, pr, n in clean:
    x = pr[0]
    print(f'  a{a:<6} gate x_{b:<6} (canonFREE={b in CANON.FREE}, val={w[b]}) pins x_{x:<6} '
          f'K={str(K)[:40]}... Kbits={abs(K).bit_length()}')

print(f'\nEMBEDDED pins (same (b,) monomial inside a big check):')
emb = [t for t in LOADS if t not in clean]
for a, b, K, pr, n in emb[:30]:
    print(f'  a{a:<6} nmon={n:<3} gate x_{b:<6} partners={pr[:8]}')

# --- who is pinned, and does it matter? ---
c2 = F2.cone([22229,22230,35758,35759,35760,35761,35762])
cc = CANON.cone([21617, 29539])
print('\n--- pinned variables vs. residual cones ---')
seen = set()
for a, b, K, pr, n in LOADS:
    for x in pr:
        if (b, x) in seen: continue
        seen.add((b, x))
        tags = []
        if x in c2: tags.append('SEVEN')
        if x in cc: tags.append('CLUSTER')
        if tags:
            print(f'  a{a:<6} gate x_{b:<6}(val {w[b]}) touches x_{x:<6} {tags} '
                  f'canonFREE={x in CANON.FREE}')

json.dump({'loads': [[a, b, str(K), pr, n] for a, b, K, pr, n in LOADS]},
          open(os.path.join(HERE,'bl_pins3.json'),'w'))
