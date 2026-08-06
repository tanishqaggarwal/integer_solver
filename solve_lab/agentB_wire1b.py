#!/usr/bin/env python3
"""Scrutinize the wire=1 failing set: categorize all 33 failing roots by gradient support,
examine the 13 unpacking equations, and redo consistency INCLUDING all failing rows as targets
(empty-gradient nonzero residual => genuinely unfixable)."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_from_solution(best)
env.tangent_linear()
res1 = env.all_root_residuals()
print(f"[b] wire=1 failing: {len(res1)}")
empty = []; withgrad = []
for i in sorted(res1):
    g = env.root_grad(i)
    (empty if not g else withgrad).append(i)
print(f"[b] failing with gradient: {len(withgrad)}: {withgrad}")
print(f"[b] failing with EMPTY gradient: {len(empty)}: {empty}")

# examine empty-gradient failing roots: what vars in their root poly, and are those free/wire/gate?
freeset = env.freeset
UNPACK = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
for i in empty[:15]:
    rp = env.root_poly[i]
    vs = set()
    for m in rp:
        vs |= set(m)
    kinds = defaultdict(int)
    for v in vs:
        if v in wire: kinds['wire'] += 1
        elif v in freeset: kinds['free'] += 1
        else: kinds['gate'] += 1
    # residual value
    print(f"  eq {i}: |terms|={len(rp)}, vars={len(vs)}, kinds={dict(kinds)}, resid_nonzero={res1[i]!=0}")

# For one empty-grad eq, print the actual monomials (var kinds) to understand structure
if empty:
    i = empty[0]
    print(f"\n[b] structure of eq {i}:")
    for m, c in list(env.root_poly[i].items())[:20]:
        tag = []
        for v in m:
            tag.append(f"x{v}({'W' if v in wire else 'F' if v in freeset else 'G'}={env.valp[v]%p if env.valp[v]%p<1000 else '...'})")
        print(f"    {c} * {'*'.join(tag) if m else '1'}")

# does the gate-cone of these depend on free inputs? check if any free input perturbation changes them
# (finite test: perturb all active free inputs a bit)
print(f"\n[b] finite test: can ANY free-input move change the empty-grad roots?")
import random; random.seed(1)
base = {v: env.valp[v] for v in freeset}
# perturb a large random subset of free inputs
hs = random.sample(sorted(freeset), 3000)
for h in hs: env.valp[h] = (env.valp[h] + random.randrange(1, p)) % p
env.forward()
changed = [i for i in empty if env.root_val(i) != res1[i]]
print(f"[b] empty-grad roots changed by random free-input move: {len(changed)}/{len(empty)}")
# restore
for v in freeset: env.valp[v] = base[v]
env.forward()
