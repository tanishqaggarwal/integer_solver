#!/usr/bin/env python3
"""POSITIVE CONTROL: validate the global mod-p solver. Perturb a constrained free input to break
some wiring eqs, then run the same global solve. It MUST heal (consistent) — proving the solver is
correct and that the CORE inconsistency is a real structural fact, not a solver artifact."""
import json, random
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}

def global_solve_step(env):
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if g: rows.append((i, g))
    pivots = {}; coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    incons = 0; piv_order = []
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0 = rows[k]; rd = dict(rd0); rhs = (-res[i]) % p if i in res else 0
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
            if rhs % p: incons += 1
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc]=({c:(v*inv)%p for c,v in rd.items()},(rhs*inv)%p); piv_order.append(pc)
    delta = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc: s = (s - v*delta.get(c,0)) % p
        delta[pc] = s % p
    return res, delta, incons

# find a constrained free input (one whose perturbation breaks wiring): from earlier, active cols.
env.set_from_solution(best); env.tangent_linear()
active = set()
for i in range(len(env.root_poly)):
    if i in env.all_root_residuals(): continue
    for c in env.root_grad(i): active.add(c)
active = sorted(active)
random.seed(11)
h = random.choice(active)
print(f"[control] perturbing constrained free input x_{h} by +12345")
env.set_from_solution(best)
env.valp[h] = (env.valp[h] + 12345) % p
env.forward()
res0 = env.all_root_residuals()
base_core = set(sorted(env.all_root_residuals()))
print(f"[control] after perturb: failing={len(res0)} (was 20 core + induced wiring breaks)")

for it in range(8):
    res, delta, incons = global_solve_step(env)
    print(f"[control] iter {it}: failing={len(res)}, incons={incons}")
    if incons:
        print(f"[control] STILL inconsistent -> the induced-break part reveals solver limits")
        # But note: the 20 CORE are also failing here; separate them
        noncore = len(res) - sum(1 for i in res if i in base_core and i in {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892})
        break
    for c, d in delta.items(): env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    if not env.all_root_residuals():
        print("[control] fully healed to 0 failing!"); break

# Cleaner control: ONLY break wiring (start from a point with NO core failing is impossible here,
# so instead: check whether the *wiring-only* subsystem (exclude 20 core rows) heals the induced break)
print("\n[control-2] wiring-only heal (exclude 20 core rows as constraints AND targets):")
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
env.set_from_solution(best)
env.valp[h] = (env.valp[h] + 12345) % p
env.forward()
for it in range(8):
    env.tangent_linear()
    res = env.all_root_residuals()
    wfail = {i: r for i, r in res.items() if i not in CORE}
    rows = []
    for i in range(len(env.root_poly)):
        if i in CORE: continue
        g = env.root_grad(i)
        if g: rows.append((i, g))
    pivots={}; coldeg=defaultdict(int)
    for _,rd in rows:
        for c in rd: coldeg[c]+=1
    incons=0; piv_order=[]
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0 = rows[k]; rd=dict(rd0); rhs=(-res[i])%p if i in res else 0
        while True:
            pc=None
            for c in rd:
                if c in pivots: pc=c; break
            if pc is None: break
            f=rd[pc]; prow,prhs=pivots[pc]
            for c,v in prow.items():
                nv=(rd.get(c,0)-f*v)%p
                if nv: rd[c]=nv
                elif c in rd: del rd[c]
            rhs=(rhs-f*prhs)%p
        if not rd:
            if rhs%p: incons+=1
            continue
        pc=min(rd,key=lambda c: coldeg.get(c,0)); inv=pow(rd[pc],p-2,p)
        pivots[pc]=({c:(v*inv)%p for c,v in rd.items()},(rhs*inv)%p); piv_order.append(pc)
    delta={}
    for pc in reversed(piv_order):
        prow,prhs=pivots[pc]; s=prhs
        for c,v in prow.items():
            if c!=pc: s=(s-v*delta.get(c,0))%p
        delta[pc]=s%p
    print(f"[control-2] iter {it}: wiring failing={len(wfail)}, incons={incons}")
    if incons: print("[control-2] wiring heal inconsistent (unexpected!)"); break
    for c,d in delta.items(): env.valp[c]=(env.valp[c]+d)%p
    env.forward()
    res2 = env.all_root_residuals()
    wfail2 = {i for i in res2 if i not in CORE}
    if not wfail2:
        print(f"[control-2] WIRING HEALED (0 wiring failures; {len(res2)} core remain) -> solver VALID"); break
