"""S12 step 19: repair the collateral while FREEZING the activation.

ac_engine.py showed the engine walks straight back to the base state -- it
repairs the collateral by undoing the activation (engine_ac_a12054_39009.json is
byte-identical to mod9118_0.json).  Forbid it from touching the activated input
and ask the real question: can the collateral be repaired with the knobs KEPT?
"""
import os, sys, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
import engine as E
P = ad.P
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
D = json.load(open(os.path.join(HERE,'ac_single.json')))
S = {int(z): r for z, r in D['res'].items()}
cands = sorted((S[z][0]['lost'], z) for z in D['pool'] if S[z][0]['knobs'] > 0)[:12]
random.seed(11)
BUDGET = int(sys.argv[1]) if len(sys.argv) > 1 else 240
print(f'freezing-repair over {len(cands)} cheapest knob activations', flush=True)
best = (B.score0, 'base', None)
for lost, z in cands:
    for tag, val in (('1', 1), ('gen', random.randrange(1 << 40, 1 << 63) | 1)):
        v = list(B.v0); v[z] = val
        A.fwd_local(v, [z])
        ch = {w for w in range(L.NVARS) if v[w] != B.v0[w]}
        sc0, newnz, newchk, l0, g0, av, nz = B.cost(v, ch)
        k0 = len(A.grad_supp(v, BAD) - supp0)
        E.FORBID = {2081, 4287, z}
        t0 = time.time()
        try:
            v2, cur = E.run(v, f'frz_{z}_{tag}', iters=40, budget=BUDGET)
        except Exception as ex:
            print(f'  x_{z} {tag}: engine failed {ex}', flush=True); continue
        av2 = L.all_atom_values(v2)
        k2 = len(A.grad_supp(v2, BAD) - supp0)
        still = v2[z] != 0
        print(f'  x_{z:<6} {tag:<4}: activated {sc0} (knobs +{k0}, {lost} eqs lost) '
              f'-> repaired {cur[0]}  activation alive {still}  knobs now +{k2}  '
              f'({time.time()-t0:.0f}s)', flush=True)
        if cur[0] > best[0]:
            best = (cur[0], f'{z}_{tag}', v2)
            T.save(v2, os.path.join(HERE, f'ac_frz_{cur[0]}.json'))
        if cur[0] > 39026:
            T.save(v2, os.path.join(HERE, 'ac_best.json'))
            print('  *** BEAT 39026 ***', flush=True)
print(f'\nBEST frozen-repair result: {best[0]} ({best[1]})')
