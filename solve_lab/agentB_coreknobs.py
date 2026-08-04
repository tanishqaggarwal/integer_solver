#!/usr/bin/env python3
"""Find free inputs that appear ONLY in the 20 core equations (pure-core knobs = safe to move,
cannot break wiring). Then build the core subsystem as a function of these knobs (with wiring
fixed at best) and attempt a nonlinear solve. Also treat the 6 override vars as candidate knobs."""
import json, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS, DEFAULT_OVERRIDE

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best); env.tangent_linear()
res0 = env.all_root_residuals(); core = set(res0)

# eqvars-based: for each free input, which equations reference it
eqvars = env.eqvars
freeset = env.freeset
free_eqs = defaultdict(set)
for i, vs in enumerate(eqvars):
    for v in vs:
        if v in freeset:
            free_eqs[v].add(i)
pure_core = [v for v, es in free_eqs.items() if es and es <= core]
print(f"[knobs] free inputs referenced in >=1 eq: {len(free_eqs)}")
print(f"[knobs] PURE-CORE free inputs (only in core eqs): {len(pure_core)}")
print(f"[knobs] pure-core list: {sorted(pure_core)}")

# which override vars are free inputs, and are they pure-core?
print("[knobs] override vars status:")
for v in DEFAULT_OVERRIDE:
    st = 'FREE' if v in freeset else 'gate'
    es = free_eqs.get(v, set())
    print(f"  x_{v}: {st}, in {len(es)} eqs, pure-core={bool(es) and es<=core}, core-eqs={sorted(es & core)}")

# quotient handles from CORE_REDUCTION
for v in [30317, 2936, 5146]:
    st = 'FREE' if v in freeset else 'gate'
    es = free_eqs.get(v, set())
    print(f"[knobs] handle x_{v}: {st}, in {len(es)} eqs, core-eqs={sorted(es & core)}, noncore={sorted(es-core)[:5]}")

# how many DISTINCT free knobs appear across core equations, and how many are 'safe'
core_free = set()
for i in core: core_free |= (eqvars[i] & freeset)
safe = [v for v in core_free if free_eqs[v] <= core]
unsafe = [v for v in core_free if not (free_eqs[v] <= core)]
print(f"[knobs] free inputs in core eqs: {len(core_free)}; safe(pure-core)={len(safe)}, "
      f"unsafe(also in wiring)={len(unsafe)}")
print(f"[knobs] unsafe core free inputs (multi-role): {sorted(unsafe)}")

# Build the core system as function of the SAFE knobs + the override vars (treat overrides as free).
# Solve mod p: vary knobs to zero the 20 core residuals. First: rank of core Jacobian wrt safe knobs.
knobs = sorted(set(safe) | set(DEFAULT_OVERRIDE))
knobs = [v for v in knobs]
print(f"[knobs] candidate solve knobs (safe pure-core + overrides): {len(knobs)}: {knobs}")

# core Jacobian wrt these knobs (treat override vars as differentiable free inputs)
# Need env where override vars are FREE. Rebuild grad with override vars added to freeset.
env.freeset = set(env.freeset) | set(DEFAULT_OVERRIDE)
env.set_from_solution(best)
# re-apply overrides values (set_from_solution used self.override for them; fine)
env.tangent_linear()
res0 = env.all_root_residuals(); coreL = sorted(res0)
J = []
for i in coreL:
    g = env.root_grad(i)
    J.append((i, {c: g[c] for c in g if c in set(knobs)}))
# rank of this core-knob Jacobian and consistency of driving residuals to 0
kidx = {k: j for j, k in enumerate(knobs)}
def rank_consistency(J, res):
    pivots = {}; incons = 0
    for i, g in sorted(J, key=lambda x: len(x[1])):
        rd = dict(g); rhs = (-res[i]) % p
        while True:
            pc = None
            for c in rd:
                if c in pivots: pc = c; break
            if pc is None: break
            f = rd[pc]; prow, prhs = pivots[pc]
            for c, v in prow.items():
                nv = (rd.get(c,0)-f*v)%p
                if nv: rd[c]=nv
                elif c in rd: del rd[c]
            rhs = (rhs-f*prhs)%p
        if not rd:
            if rhs % p: incons += 1
            continue
        pc = min(rd); inv = pow(rd[pc],p-2,p)
        pivots[pc]=({c:(v*inv)%p for c,v in rd.items()},(rhs*inv)%p)
    return len(pivots), incons
rk, inc = rank_consistency(J, res0)
print(f"[knobs] core Jacobian wrt {len(knobs)} knobs: rank={rk}, inconsistent={inc} / {len(coreL)} core eqs")
print(f"[knobs] => core linearly {'SOLVABLE' if inc==0 else 'UNREACHABLE'} with these knobs")
