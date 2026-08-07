"""Rank knobs by (dA,dB) direction and by the equation-cost of the atoms they break."""
import json, sys, time, collections
import dlib as L
import engine2 as E
import adv3
P = L.P
A_, B_ = 35389, 6671

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json'))
A0, B0 = st0.v[A_] % P, st0.v[B_] % P
base_nz = set(st0.nz())
print('base', st0.score, sorted(base_nz))

rows = []
t0 = time.time()
for k, u in enumerate(sorted(L.freeset)):
    st = st0.clone()
    st.apply({u: st.v[u] + 1})
    adv3.sweep(st, rounds=6)
    a, b = st.v[A_] % P, st.v[B_] % P
    da, db = (a - A0) % P, (b - B0) % P
    if da == 0 and db == 0:
        continue
    nz = set(st.nz())
    new = sorted(nz - base_nz)
    eqs = set()
    for c in new:
        eqs |= set(L.atom2eq.get(c, {}))
    rows.append((u, da, db, new, len(eqs), st.score))
    if k % 2000 == 0:
        print('  ', k, f'{time.time()-t0:.0f}s', len(rows), flush=True)

print('movers of (A,B):', len(rows))
rows.sort(key=lambda r: r[4])
print('cheapest by #equations of newly-broken atoms:')
for u, da, db, new, ne, sc in rows[:40]:
    print(f'  x_{u:<7} neweqs={ne:<4} score={sc} brokenatoms={new}')
json.dump([[u, str(da), str(db), new, ne, sc] for u, da, db, new, ne, sc in rows],
          open('scanpairs.json', 'w'))

# --- direction analysis: group by (dA:dB) ratio mod p ---
def ratio(da, db):
    if db == 0:
        return ('inf',)
    return (da * pow(db, P - 2, P) % P,)

byratio = collections.defaultdict(list)
for u, da, db, new, ne, sc in rows:
    byratio[ratio(da, db)].append((ne, u))
print('distinct (dA:dB) directions:', len(byratio))
best = sorted(((min(v), k) for k, v in byratio.items()))[:15]
for (ne, u), k in best:
    print(f'   direction {str(k)[:40]:<42} cheapest x_{u} neweqs={ne}')
