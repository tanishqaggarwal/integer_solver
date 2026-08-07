"""S12 step 13: exhaustive pair sweep AGAIN with generic (random 64-bit) values,
to rule out that value 1 caused accidental cancellations, and with the correct
control: pair knobs vs the better of the two SINGLES at the same values."""
import os, sys, json, time, itertools, collections, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
D = json.load(open(os.path.join(HERE,'ac_single.json')))
pool = D['pool']
random.seed(2718)
VAL = {z: random.randrange(1 << 40, 1 << 63) | 1 for z in pool}
single = {}
for z in pool:
    v = list(B.v0); v[z] = VAL[z]
    A.fwd_local(v, [z])
    single[z] = len(A.grad_supp(v, BAD) - supp0)
print(f'singles at generic values: knobs>0 for {sum(1 for z in pool if single[z])} of {len(pool)}; '
      f'max {max(single.values())}', flush=True)
t0 = time.time(); syn = []; hits = 0; mx = 0
pairs = list(itertools.combinations(pool, 2))
for i, (u, z) in enumerate(pairs):
    v = list(B.v0); v[u] = VAL[u]; v[z] = VAL[z]
    A.fwd_local(v, [u, z])
    k = len(A.grad_supp(v, BAD) - supp0)
    if k: hits += 1
    mx = max(mx, k)
    m = max(single[u], single[z])
    if k > m:
        changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
        sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
        syn.append((k - m, k, m, len(lost), u, z, sc))
    if i % 4000 == 0: print(f'  {i}/{len(pairs)} hits {hits} max {mx} syn {len(syn)} ({time.time()-t0:.0f}s)', flush=True)
print(f'\nswept {len(pairs)} generic-value pairs in {time.time()-t0:.0f}s')
print(f'  pairs growing the support: {hits};  max knobs from any pair: {mx}')
print(f'  GENUINELY SECOND-ORDER pairs (beat both singles): {len(syn)}')
syn.sort(key=lambda t: (-t[0], t[3]))
for t in syn[:25]: print(f'   extra={t[0]} knobs={t[1]} best_single={t[2]} eqs_lost={t[3]} x_{t[4]},x_{t[5]} score={t[6]}')
json.dump({'syn': syn, 'max': mx, 'hits': hits}, open(os.path.join(HERE,'ac_pairs2.json'),'w'))
