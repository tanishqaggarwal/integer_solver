#!/usr/bin/env python3
"""Test the CORE_REDUCTION cancellation conditions and whether the override knobs help.
Checks x_24468-x_13682, x_18956-x_37892 (C1/C2 cancels), x_29322, x_3558 residues at best; then
tries setting each override to its 'cancellation' target and moving x_30213 (the safe wiring-null
knob), counting failing equations after re-forward."""
import json
from agentB_setup import load, Env, p
data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
def V(v): return env.valp[v]
def nfail(): return len(env.all_root_residuals())
print(f"[cancel] best failing = {nfail()}")
print(f"[cancel] x_24468 - x_13682 = {(V(24468)-V(13682))%p}")
print(f"[cancel] x_18956 - x_37892 = {(V(18956)-V(37892))%p}")
print(f"[cancel] x_29322 = x_14853 - x_12186 = {(V(14853)-V(12186))%p}")
print(f"[cancel] x_3558  = x_24908 - x_16742 = {(V(24908)-V(16742))%p}")
print(f"[cancel] base gate S=x_35389={V(35389)}, T=x_6671={V(6671)}")

# try setting overrides to cancellation targets
for name, ov, tgt in [("x_24468:=x_13682", 24468, 13682), ("x_18956:=x_37892", 18956, 37892)]:
    env.set_from_solution(best)
    env.valp[ov] = V(tgt)
    env.forward()
    print(f"[cancel] set {name}: failing={nfail()}")

# both together
env.set_from_solution(best)
env.valp[24468] = env.valp[13682]; env.valp[18956] = env.valp[37892]
env.forward()
print(f"[cancel] set BOTH cancels: failing={nfail()}")

# move x_30213 (safe knob) across values; count failing (wiring stays satisfied)
env.set_from_solution(best)
base30213 = env.valp[30213]
import random; random.seed(0)
print("[cancel] sweeping x_30213 (safe wiring-null knob):")
best_seen = (20, None)
for tag, val in [("0",0),("1",1),("C1",best.get(22162,env.valp[22162])),
                 ("rand1",random.randrange(p)),("rand2",random.randrange(p)),
                 ("x_13682", None)]:
    env.set_from_solution(best)
    v = env.valp[13682] if tag=="x_13682" else val
    env.valp[30213] = v % p; env.forward()
    f = nfail()
    if f < best_seen[0]: best_seen = (f, tag)
    print(f"   x_30213={tag}: failing={f}")
print(f"[cancel] best over sweep: {best_seen}")
env.set_from_solution(best)
