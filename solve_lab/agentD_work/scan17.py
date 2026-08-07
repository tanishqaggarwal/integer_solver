"""At D_39017 the only nonzero atoms are a688, a1618, a40608.  Scan every free
input for moves that keep the score and change those three values -> the lattice
of reachable residuals, then price every subset of the 16 failing equations."""
import json, sys, time, itertools, collections
import dlib as L
import engine2 as E
import adv3
P = L.P
RES = [688, 1618, 40608]

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_39017.json'))
print('base', st0.score, st0.nz())
base = [st0.av[a] for a in RES]
fails = sorted(st0.fail)
print('failing eqs', fails)
for i in fails:
    m, sq, co = L.eq_atoms[i]
    print(f'  eq{i}: mult={m} sq={sq} atoms={ {a: c for a, c in co.items()} }')

good = []
t0 = time.time()
for k, u in enumerate(sorted(L.freeset)):
    st = st0.clone()
    st.apply({u: st.v[u] + 1})
    d = [st.av[a] - base[i] for i, a in enumerate(RES)]
    if st.score >= st0.score and any(d):
        good.append((u, st.score, d))
    if k % 2000 == 0:
        print('   ', k, f'{time.time()-t0:.0f}s', len(good), flush=True)
print('cost-free residual movers:', len(good))
for u, s, d in good[:40]:
    print('   x_%d score=%d delta=%s' % (u, s, [x for x in d]))
json.dump([[u, s, [str(x) for x in d]] for u, s, d in good], open('scan17.json', 'w'))
