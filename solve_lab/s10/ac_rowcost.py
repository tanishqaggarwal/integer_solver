"""S12 step 16: the decisive ratio -- do the new COLUMNS an activation buys add
fewer ROWS than columns?

For every activation that grows the cluster's gradient support, compute each new
knob's jac_column and count how many of the checks it touches are NOT already
rows of the base closure.  columns_gained / rows_gained > 1 is the only regime in
which a kernel can ever open.
"""
import os, sys, json, time, random, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
from fwdad import jac_column
P = ad.P; FORBID = {2081, 4287}
CHECKS = sorted(a for a in range(L.NA) if a not in L.atom_out)
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
S = json.load(open(os.path.join(HERE,'ac_sacr.json')))
BASEROWS = set(S['rows']); BASECOLS = set(S['Us'])
print(f'base closure {len(BASEROWS)} rows x {len(BASECOLS)} cols', flush=True)
D = json.load(open(os.path.join(HERE,'ac_single.json')))
pool = D['pool']; SG = {int(z): r for z, r in D['res'].items()}
cand = [z for z in pool if SG[z][0]['knobs'] > 0 or SG[z][2]['knobs'] > 0]
print(f'activations that grow the support: {len(cand)}', flush=True)
random.seed(99)
out = []
t0 = time.time()
for i, z in enumerate(cand):
    for tagv, val in (('1', 1), ('gen', random.randrange(1 << 40, 1 << 63) | 1)):
        v = list(B.v0); v[z] = val
        A.fwd_local(v, [z])
        vm = [x % P for x in v]
        supp = A.grad_supp(v, BAD)
        newk = sorted(supp - supp0); lostk = sorted(supp0 - supp)
        if not newk: continue
        changed = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
        sc, newnz, newchk, lost, gained, av, nz = B.cost(v, changed)
        newrows = set()
        for u in newk:
            newrows |= set(jac_column(u, v, vm, CHECKS))
        nr = newrows - BASEROWS
        nc = len(newk)
        out.append((nc/(len(nr)+1e-9) if nr else float('inf'), nc, len(nr), len(lostk),
                    len(lost), z, tagv, sc, sorted(nr)[:10]))
    if i % 15 == 0: print(f'  {i}/{len(cand)} ({time.time()-t0:.0f}s)', flush=True)
out.sort(reverse=True)
print(f'\nranked by columns-gained per NEW row ({time.time()-t0:.0f}s):')
print('   ratio  newcols  newrows  cols_lost  eqs_lost   x      val   score')
for r in out[:30]:
    print(f'  {r[0]:>7.3f}  {r[1]:>6}  {r[2]:>7}  {r[3]:>8}  {r[4]:>8}  x_{r[5]:<6} {r[6]:<4} {r[7]}')
good = [r for r in out if r[2] == 0]
print(f'\nactivations whose new knobs add NO new closure row: {len(good)}')
for r in good[:20]: print(f'  {r}')
json.dump([[r[0] if r[0] != float("inf") else -1] + list(r[1:]) for r in out],
          open(os.path.join(HERE,'ac_rowcost.json'),'w'))
