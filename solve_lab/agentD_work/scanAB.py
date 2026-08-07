"""Scan every free input: perturb, re-solve the advice congruences, record
(score, A mod p, B mod p).  A = x_35389, B = x_6671.  At D_adv the whole instance
is satisfied iff A == B == 0 (mod p)."""
import json, sys, time, collections, random
import dlib as L
import engine2 as E
import adv3
P = L.P
A_, B_ = 35389, 6671

base = L.load(sys.argv[1] if len(sys.argv) > 1 else 'D_adv.json')
DELTAS = [int(x) for x in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['1'])]
st0 = E.St(base)
print('base score', st0.score, 'nz', st0.nz())
print('A =', st0.v[A_] % P)
print('B =', st0.v[B_] % P)
A0, B0 = st0.v[A_] % P, st0.v[B_] % P

out = []
t0 = time.time()
frees = sorted(L.freeset)
for k, u in enumerate(frees):
    for d in DELTAS:
        st = st0.clone()
        st.apply({u: st.v[u] + d})
        adv3.sweep(st, rounds=6)
        a, b = st.v[A_] % P, st.v[B_] % P
        if st.score != st0.score or a != A0 or b != B0:
            out.append((u, d, st.score, a != A0, b != B0))
    if k % 500 == 0:
        print(f'  {k}/{len(frees)}  {time.time()-t0:.0f}s  interesting={len(out)}', flush=True)
print(f'scan done {time.time()-t0:.0f}s')
neutral = [r for r in out if r[2] >= st0.score and (r[3] or r[4])]
print('MOVES THAT CHANGE (A,B) AT NO COST:', len(neutral))
for r in neutral[:50]:
    print('   ', r)
cheap = sorted([r for r in out if (r[3] or r[4])], key=lambda r: -r[2])
print('top movers of (A,B) by score:')
for r in cheap[:40]:
    print('   ', r)
json.dump(out, open('scanAB.json', 'w'))
