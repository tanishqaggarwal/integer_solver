#!/usr/bin/env python3
"""Understand the 3 remaining core eqs after Route A. Are they M1/M2/M3? Do they have
private handles? Check dependence on the changed slack x_16742/x_14853."""
import json, re
from agentA_harness import (p, load_solution, forward, eval_fails, NEQ, lines, eqvars,
                            freeinp, backward_cone, gates, definer, order)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)

def eqval(v, c):
    code = compile(re.sub(r'x_(\d+)', r'v[\1]', lines[c].rsplit('=', 1)[0]), '<e>', 'eval')
    return eval(code, {'__builtins__': {}, 'v': v})

# baseline values of the 3 core
core3 = [11854, 29437, 32916]
print("=== baseline values of remaining-3 core ===")
for c in core3:
    print(f" eq {c}: {eqval(base,c)}")

# Route A
v = base[:]; v[16742] = base[24908]; v[14853] = base[12186]; forward(v)
for hv, expr in [(30317, -(v[11150])//p if v[11150] % p == 0 else 0),
                 (2936, (537773*v[37758])//p if (537773*v[37758]) % p == 0 else 0),
                 (5146, v[25739]//(6672769*p) if v[25739] % (6672769*p) == 0 else 0)]:
    v[hv] = expr
print("\n=== Route A values of remaining-3 core ===")
for c in core3:
    val = eqval(v, c)
    print(f" eq {c}: {val}   mod p={val%p}")

# which free inputs (private-ish) does each core3 depend on, that are NOT in other equations?
# count how many equations each free input in core3's cone appears in
gdef_vids = {t: gates[definer[t]][2] for t in order}
def free_cone(vs):
    seen=set(); st=list(vs)
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gdef_vids.get(u,()):
            if w not in seen: st.append(w)
    return set(u for u in seen if u in freeinp)
var_eqcount = {}
from collections import Counter
allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x]+=1

print("\n=== core3 private free-input handles (appear in few eqs) ===")
for c in core3:
    fc = free_cone(eqvars[c])
    # direct free vars in the eq
    directfree = [x for x in eqvars[c] if x in freeinp]
    rare = [(x, allvarcount[x]) for x in directfree if allvarcount[x] <= 3]
    print(f" eq {c}: direct free vars={sorted(directfree)}")
    print(f"    rare (<=3 eqs): {rare}")
    print(f"    depends on x_16742? {16742 in fc}  x_14853? {14853 in fc}  x_12186? {12186 in fc}")
