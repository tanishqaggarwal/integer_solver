#!/usr/bin/env python3
"""Analyze the 23 broken + 3 remaining-core after Route A. Find healing free inputs."""
import json
from collections import defaultdict, Counter
from agentA_harness import (p, load_solution, forward, eval_fails, NEQ, order,
                            gates, definer, freeinp, eqvars, lines, backward_cone)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
F0 = set(eval_fails(base))

# critical vars we must NOT change (they define S,T=0)
crit = {16742, 14853, 12186, 24908, 3558, 29322}
# add S,T sub-cones (a=x_33469, etc) - free inputs feeding them
for r in [33469, 29322, 27713, 1326, 3558, 24908]:
    _, fr = backward_cone(r); crit |= fr

v = base[:]
v[16742] = base[24908]; v[14853] = base[12186]
forward(v)
for hv, expr in [(30317, -(v[11150])//p if v[11150] % p == 0 else 0),
                 (2936, (537773*v[37758])//p if (537773*v[37758]) % p == 0 else 0),
                 (5146, v[25739]//(6672769*p) if v[25739] % (6672769*p) == 0 else 0)]:
    v[hv] = expr
F = set(eval_fails(v))
broke = sorted(F - F0)
remcore = sorted(F & F0)   # still-failing originals
print(f"broke {len(broke)}: {broke}")
print(f"remaining core {len(remcore)}: {remcore}")

# For remaining core, show equations
import re as _re
def eqval(c):
    code = compile(_re.sub(r'x_(\d+)', r'v[\1]', lines[c].rsplit('=', 1)[0]), '<e>', 'eval')
    return eval(code, {'__builtins__': {}, 'v': v})
print("\n=== remaining core eqs ===")
for c in remcore:
    print(f" eq {c}: nvars={len(eqvars[c])} val={eqval(c)}")
    print(f"   vars: {sorted(eqvars[c])}")

# For broken, find free-input cones and healing candidates (free inputs feeding the eq, not critical)
gdef_vids = {t: gates[definer[t]][2] for t in order}
def free_cone(varset):
    seen=set(); st=list(varset)
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gdef_vids.get(u,()):
            if w not in seen: st.append(w)
    return set(u for u in seen if u in freeinp)

print("\n=== broken eqs: healing candidates ===")
heal_knobs = Counter()
brk_info = {}
for c in broke:
    fc = free_cone(eqvars[c])
    cand = fc - crit
    brk_info[c] = sorted(cand)
    for k in cand: heal_knobs[k]+=1
    print(f" eq {c}: {len(fc)} free-deps, {len(cand)} non-critical healing candidates")
print(f"\nmost common healing knobs (free input -> #broken eqs it feeds):")
for k,ct in heal_knobs.most_common(30):
    print(f"   x_{k}: {ct}  (base val {base[k]})")
json.dump({'broke':broke,'remcore':remcore,'brk_info':brk_info,
           'heal_knobs':dict(heal_knobs)}, open('agentA_heal_analyze.json','w'))
