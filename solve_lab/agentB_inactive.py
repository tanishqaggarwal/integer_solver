#!/usr/bin/env python3
"""Decisive nonlinearity test. At best: split free inputs into ACTIVE (nonzero root-gradient)
and INACTIVE. Perturb ONLY inactive columns by finite random amounts; forward-eval; measure
(a) how many wiring eqs break, (b) whether the 20 core residuals change. If inactive moves keep
wiring satisfied but move the core -> a nonlinear handle exists. Also test: do the specific
control free inputs (x_14853,x_16742,x_12186,x_24908) reach the core, and at what cost."""
import json, random, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best); env.tangent_linear()
res0 = env.all_root_residuals(); core = sorted(res0)

active = set()
for i in range(len(env.root_poly)):
    for c in env.root_grad(i): active.add(c)
inactive = [v for v in env.freeinp if v not in active]
print(f"[inact] active cols={len(active)}, inactive free inputs={len(inactive)}")

# baseline residual vector over all eqs
def resid_vec():
    return {i: env.root_val(i) for i in range(len(env.root_poly)) if env.root_val(i)}

base = {v: env.valp[v] for v in env.freeset}
random.seed(7)
for trial, npert in enumerate([100, 1000, len(inactive)]):
    for v in env.freeset: env.valp[v] = base[v]
    hs = random.sample(inactive, min(npert, len(inactive)))
    for h in hs: env.valp[h] = random.randrange(p)
    env.forward()
    rv = resid_vec()
    fails = set(rv)
    new_break = fails - set(core)
    core_changed = sum(1 for i in core if rv.get(i, 0) != res0[i])
    print(f"[inact] perturb {len(hs)} inactive cols: total fail={len(fails)}, "
          f"NEW wiring breaks={len(new_break)}, core residuals changed={core_changed}/20")
# restore
for v in env.freeset: env.valp[v] = base[v]
env.forward()

# targeted: do the control differences reach S,T? Try the CORE_REDUCTION finite fix.
# x_29322 = x_14853 - x_12186 ; x_3558 = x_24908 - x_16742
def valof(v): return env.valp[v]
print("\n[inact] control gate residues at best:")
for name, expr in [("x_29322", (14853, 12186)), ("x_3558", (24908, 16742))]:
    a, b = expr
    print(f"  {name} = x_{a}-x_{b} = {(valof(a)-valof(b))%p}  (a in free={a in env.freeset}, b in free={b in env.freeset})")
# apply finite fix and count breakage
x29322 = (valof(14853) - valof(12186)) % p
x3558 = (valof(24908) - valof(16742)) % p
env.valp[14853] = (env.valp[14853] - x29322) % p
env.valp[16742] = (env.valp[16742] + x3558) % p
env.forward()
rv = resid_vec()
print(f"[inact] after finite fix (x_14853,x_16742): total fail(modp)={len(rv)}, "
      f"core still failing={sum(1 for i in core if i in rv)}, new breaks={len(set(rv)-set(core))}")
