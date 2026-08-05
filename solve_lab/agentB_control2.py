#!/usr/bin/env python3
"""Fast positive control (forward-elimination consistency only).
(A) Perturb a constrained free input -> breaks wiring. The WIRING-ONLY global system (exclude the
    20 core) MUST be consistent (the undo delta heals it). Validates solver + confirms the CORE is
    the sole obstruction.
(B) Same test but INCLUDE the perturbed break as a target together with the 20 core: shows the 20
    core stay inconsistent regardless."""
import json, random
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}

def consistency(env, exclude_rows, extra_target_ok=True):
    """Forward-eliminate all nonzero rows (excluding `exclude_rows` from being constraints);
    RHS = -residual. Return (rank, incons_count, incons_list)."""
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        if i in exclude_rows: continue
        g = env.root_grad(i)
        if g: rows.append((i, g, (-res[i]) % p if i in res else 0))
    pivots = {}; coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    incons = 0; ilist = []
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0, rhs = rows[k]; rd = dict(rd0)
        while True:
            pc = None
            for c in rd:
                if c in pivots: pc = c; break
            if pc is None: break
            f = rd[pc]; prow, prhs = pivots[pc]
            for c, v in prow.items():
                nv = (rd.get(c,0)-f*v) % p
                if nv: rd[c]=nv
                elif c in rd: del rd[c]
            rhs = (rhs-f*prhs) % p
        if not rd:
            if rhs % p: incons += 1; ilist.append(i)
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p)
    return len(pivots), incons, ilist

# pick a constrained free input
env.set_from_solution(best); env.tangent_linear()
active = set()
for i in range(len(env.root_poly)):
    if i in CORE: continue
    for c in env.root_grad(i): active.add(c)
random.seed(11); h = random.choice(sorted(active))

# baseline (no perturbation): wiring-only excluding core must be consistent (it's already satisfied)
env.set_from_solution(best)
rk, inc, il = consistency(env, exclude_rows=CORE)
print(f"[ctrl-A0] baseline, exclude core: rank={rk}, incons={inc}  (expect 0)")

# perturb h -> break wiring; wiring-only heal must be consistent
env.set_from_solution(best); env.valp[h] = (env.valp[h] + 12345) % p; env.forward()
res = env.all_root_residuals()
print(f"[ctrl-A1] perturbed x_{h}: total failing={len(res)} ({len(res)-20} induced wiring breaks + 20 core)")
rk, inc, il = consistency(env, exclude_rows=CORE)
print(f"[ctrl-A1] wiring-only heal (exclude core): rank={rk}, incons={inc}  -> "
      f"{'CONSISTENT: solver heals induced break (VALID)' if inc==0 else 'INCONSISTENT'}")

# include core as targets too: core stays inconsistent
env.set_from_solution(best)
rk, inc, il = consistency(env, exclude_rows=set())
print(f"[ctrl-B] full system incl core (best): rank={rk}, incons={inc}, incons_rows={il[:8]} "
      f"-> core is the obstruction")
