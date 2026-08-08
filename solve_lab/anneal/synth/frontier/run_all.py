#!/usr/bin/env python3
"""run_all.py -- consolidated frontier measurement, priority-ordered, each section
checkpointed to its own JSON so partial runs are still usable."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import frontier as F

T0 = time.time()
def stamp(msg): print(f"\n##### {msg}  (+{time.time()-T0:.0f}s) #####", flush=True)


# ---- Section A: CORE solver x s x effort table (baseline wallace) ----------
stamp("SECTION A: core solver x s x effort")
F.sweep(sizes=[4, 5, 6, 7, 8],
        solvers=['sa', 'pt', 'tabu', 'sb'],
        efforts=['low', 'mid', 'high', 'extreme'],
        encoder='baseline', mode='wallace',
        n_caps=dict(sa=25, pt=12, tabu=250, sb=300),
        time_cap=15.0, out_json=os.path.join(HERE, 'core_baseline.json'), label='core')

# ---- Section B: combined solvers (can they break the wall?) ----------------
stamp("SECTION B: combos sb_tabu, pt_wide at the wall")
F.sweep(sizes=[6, 7, 8],
        solvers=['sb_tabu', 'pt_wide'],
        efforts=['high', 'extreme'],
        encoder='baseline', mode='wallace',
        n_caps=dict(sb_tabu=200, pt_wide=8),
        time_cap=20.0, out_json=os.path.join(HERE, 'combos.json'), label='combo')

# ---- Section C: clamp one operand (a known) -- the "one operand known" F ----
stamp("SECTION C: clamp operand a (one operand known)")
F.sweep(sizes=[4, 5, 6, 7, 8, 10, 12],
        solvers=['tabu', 'pt'],
        efforts=['mid', 'high'],
        encoder='baseline', mode='wallace', clamp_which='a',
        n_caps=dict(tabu=250, pt=12),
        time_cap=15.0, out_json=os.path.join(HERE, 'clamp_a.json'), label='clampA')

# ---- Section D: squeeze encoder vs baseline, head to head ------------------
stamp("SECTION D: squeeze (karatsuba+NAF) encoder")
for mode, mult in [('wallace', 'karatsuba'), ('dadda', 'karatsuba')]:
    F.sweep(sizes=[4, 5, 6, 7, 8],
            solvers=['tabu', 'sb'],
            efforts=['mid', 'high'],
            encoder='squeeze',
            squeeze_kw=dict(mult=mult, red='naf', leaf=8, mode=mode),
            n_caps=dict(tabu=250, sb=300),
            time_cap=15.0,
            out_json=os.path.join(HERE, f'squeeze_{mode}.json'),
            label=f'sq-{mode}')

# ---- Section E: penalty-weight / AND-weight (W_and) + binary mode ----------
stamp("SECTION E: W_and tuning + binary mode")
res = {}
for mode, W_and in [('wallace', None), ('wallace', 2), ('wallace', 8),
                    ('wallace', 32), ('wallace', 128), ('binary', None)]:
    r = F.sweep(sizes=[5, 6, 7], solvers=['tabu'], efforts=['mid', 'high'],
                encoder='baseline', mode=mode, W_and=W_and,
                n_caps=dict(tabu=250), time_cap=15.0,
                label=f'{mode}/Wand={W_and}')
    for k, v in r.items():
        v['mode'] = mode; v['W_and'] = W_and
        res[f"{mode}|{W_and}|{k}"] = v
    json.dump(res, open(os.path.join(HERE, 'wand.json'), 'w'), indent=1)

stamp("ALL DONE")
