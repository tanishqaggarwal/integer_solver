#!/usr/bin/env python3
"""run_rest.py -- sections not finished by run_all2 (killed in Section B when
pt_wide/extreme ran 300s+ per restart).  Fast tabu/sb only; pt_wide hard-capped."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import frontier as F

T0 = time.time()
def stamp(m): print(f"\n##### {m}  (+{time.time()-T0:.0f}s) #####", flush=True)
def save(d, name): json.dump(d, open(os.path.join(HERE, name), 'w'), indent=1)
def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else {}

# ---- B (finish): combos, fast sb_tabu + hard-capped pt_wide ----------------
stamp("SECTION B finish: combos")
combos = load('combos.json')
combos.update(F.sweep(sizes=[6, 7, 8], solvers=['sb_tabu'],
                      efforts=['high', 'extreme'], encoder='baseline', mode='wallace',
                      n_caps=dict(sb_tabu=150), time_cap=12.0, label='combo'))
save(combos, 'combos.json')
# pt_wide is ~100-330s/run: allow at most 2 runs at s=7 high as a spot-check
combos.update(F.sweep(sizes=[7], solvers=['pt_wide'], efforts=['high'],
                      encoder='baseline', mode='wallace',
                      n_caps=dict(pt_wide=2), time_cap=90.0, label='combo-ptw'))
save(combos, 'combos.json')

# ---- C: clamp one operand (one operand known) ------------------------------
stamp("SECTION C: clamp a")
clamp = {}
clamp.update(F.sweep(sizes=[4, 5, 6, 7, 8, 10, 12], solvers=['tabu'],
                     efforts=['mid', 'high', 'extreme'], encoder='baseline',
                     mode='wallace', clamp_which='a',
                     n_caps=dict(tabu=250), time_cap=8.0, label='clampA'))
save(clamp, 'clamp_a.json')

# ---- D: squeeze encoder vs baseline ----------------------------------------
stamp("SECTION D: squeeze")
for mode, mult in [('wallace', 'karatsuba'), ('dadda', 'karatsuba')]:
    F.sweep(sizes=[4, 5, 6, 7, 8], solvers=['tabu', 'sb'],
            efforts=['mid', 'high'], encoder='squeeze',
            squeeze_kw=dict(mult=mult, red='naf', leaf=8, mode=mode),
            n_caps=dict(tabu=250, sb=300), time_cap=8.0,
            out_json=os.path.join(HERE, f'squeeze_{mode}.json'), label=f'sq-{mode}')

# ---- E: W_and tuning + binary mode -----------------------------------------
stamp("SECTION E: W_and + binary")
res = {}
for mode, W_and in [('wallace', None), ('wallace', 2), ('wallace', 8),
                    ('wallace', 32), ('wallace', 128), ('binary', None)]:
    r = F.sweep(sizes=[5, 6, 7], solvers=['tabu'], efforts=['mid', 'high'],
                encoder='baseline', mode=mode, W_and=W_and,
                n_caps=dict(tabu=250), time_cap=8.0, label=f'{mode}/W={W_and}')
    for k, v in r.items():
        v['mode'] = mode; v['W_and'] = W_and
        res[f"{mode}|{W_and}|{k}"] = v
    save(res, 'wand.json')

stamp("ALL DONE")
