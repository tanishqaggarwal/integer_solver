#!/usr/bin/env python3
"""run_all2.py -- leaner, priority-ordered frontier measurement.
Fast solvers (tabu, sb) carry the dense effort sweep; slow sa/pt run only where a
rate is affordable, with a couple of high/extreme spot-checks to show big budgets
still fail.  Each section checkpoints its own JSON."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import frontier as F

T0 = time.time()
def stamp(m): print(f"\n##### {m}  (+{time.time()-T0:.0f}s) #####", flush=True)
def save(d, name): json.dump(d, open(os.path.join(HERE, name), 'w'), indent=1)


# ============ Section A: core solver x s x effort (baseline wallace) ========
stamp("SECTION A: core")
core = {}
# fast solvers: full effort sweep, s=4..9
core.update(F.sweep(sizes=[4, 5, 6, 7, 8, 9], solvers=['tabu', 'sb'],
                    efforts=['low', 'mid', 'high', 'extreme'],
                    encoder='baseline', mode='wallace',
                    n_caps=dict(tabu=250, sb=300), time_cap=8.0, label='core'))
save(core, 'core_baseline.json')
# slow solvers: affordable efforts (low, mid) across the frontier
core.update(F.sweep(sizes=[4, 5, 6, 7, 8, 9], solvers=['sa', 'pt'],
                    efforts=['low', 'mid'],
                    encoder='baseline', mode='wallace',
                    n_caps=dict(sa=20, pt=12), time_cap=12.0, label='core-slow'))
save(core, 'core_baseline.json')
# slow solvers: high/extreme spot-check at the wall (few runs, best-E evidence)
core.update(F.sweep(sizes=[6, 7], solvers=['sa', 'pt'],
                    efforts=['high', 'extreme'],
                    encoder='baseline', mode='wallace',
                    n_caps=dict(sa=6, pt=4), time_cap=18.0, label='core-slow-hi'))
save(core, 'core_baseline.json')

# ============ Section B: combined solvers -- break the wall? =================
stamp("SECTION B: combos")
F.sweep(sizes=[6, 7, 8], solvers=['sb_tabu', 'pt_wide'],
        efforts=['high', 'extreme'], encoder='baseline', mode='wallace',
        n_caps=dict(sb_tabu=200, pt_wide=6), time_cap=12.0,
        out_json=os.path.join(HERE, 'combos.json'), label='combo')

# ============ Section C: clamp one operand (one operand known) ==============
stamp("SECTION C: clamp a")
clamp = {}
clamp.update(F.sweep(sizes=[4, 5, 6, 7, 8, 10, 12], solvers=['tabu'],
                     efforts=['mid', 'high', 'extreme'], encoder='baseline',
                     mode='wallace', clamp_which='a',
                     n_caps=dict(tabu=250), time_cap=8.0, label='clampA'))
save(clamp, 'clamp_a.json')
clamp.update(F.sweep(sizes=[6, 7, 8], solvers=['pt'], efforts=['mid'],
                     encoder='baseline', mode='wallace', clamp_which='a',
                     n_caps=dict(pt=10), time_cap=12.0, label='clampA-pt'))
save(clamp, 'clamp_a.json')

# ============ Section D: squeeze encoder vs baseline ========================
stamp("SECTION D: squeeze")
for mode, mult in [('wallace', 'karatsuba'), ('dadda', 'karatsuba')]:
    F.sweep(sizes=[4, 5, 6, 7, 8], solvers=['tabu', 'sb'],
            efforts=['mid', 'high'], encoder='squeeze',
            squeeze_kw=dict(mult=mult, red='naf', leaf=8, mode=mode),
            n_caps=dict(tabu=250, sb=300), time_cap=8.0,
            out_json=os.path.join(HERE, f'squeeze_{mode}.json'), label=f'sq-{mode}')

# ============ Section E: W_and tuning + binary mode ========================
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
