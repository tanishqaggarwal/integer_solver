#!/usr/bin/env python3
import json, time, sys
from agentB_setup import load, Env, p, NVARS

t0 = time.time()
data = load()
env = Env(data)
print(f"free inputs: {len(env.freeinp)}, gates(order): {len(env.order)}, eqs: {len(env.root_poly)}", file=sys.stderr)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
res = env.all_root_residuals()
print(f"[probe] nonzero roots mod p: {len(res)}  in {time.time()-t0:.1f}s")
print("indices:", sorted(res)[:40])

# measure gradient support growth via tangent-linear pass wrt ALL free inputs
# grad[v] = dict{freeinput: coef}. free input j: {j:1}. Track support sizes.
vp = env.valp
grad = [None] * NVARS
for v in env.freeset:
    grad[v] = {v: 1}
maxsup = 0; totsup = 0; ng = 0
t1 = time.time()
for t, pol in env.gate_poly:
    g = defaultdict = {}
    for m, c in pol.items():
        if len(m) == 0:
            continue
        elif len(m) == 1:
            v = m[0]
            gv = grad[v]
            if gv:
                for k, cc in gv.items():
                    g[k] = (g.get(k, 0) + c * cc) % p
        elif len(m) == 2:
            a, b = m
            if a == b:
                # c * a^2 -> 2 c a * grad[a]
                fac = (2 * c * vp[a]) % p
                gv = grad[a]
                if gv and fac:
                    for k, cc in gv.items():
                        g[k] = (g.get(k, 0) + fac * cc) % p
            else:
                fa = (c * vp[b]) % p; fb = (c * vp[a]) % p
                ga = grad[a]; gb = grad[b]
                if ga and fa:
                    for k, cc in ga.items():
                        g[k] = (g.get(k, 0) + fa * cc) % p
                if gb and fb:
                    for k, cc in gb.items():
                        g[k] = (g.get(k, 0) + fb * cc) % p
    g = {k: v2 for k, v2 in g.items() if v2}
    grad[t] = g
    s = len(g)
    if s > maxsup: maxsup = s
    totsup += s; ng += 1
print(f"[probe] tangent-linear pass over {ng} gates in {time.time()-t1:.1f}s")
print(f"[probe] gradient support: max={maxsup}, mean={totsup/max(ng,1):.1f}")

# root gradient support for the nonzero (core) roots
from collections import defaultdict as dd
def root_grad(i):
    g = {}
    for m, c in env.root_poly[i].items():
        if len(m) == 0: continue
        elif len(m) == 1:
            v = m[0]; gv = grad[v]
            if gv:
                for k, cc in gv.items(): g[k] = (g.get(k, 0) + c * cc) % p
        elif len(m) == 2:
            a, b = m
            if a == b:
                fac = (2 * c * vp[a]) % p; gv = grad[a]
                if gv and fac:
                    for k, cc in gv.items(): g[k] = (g.get(k, 0) + fac * cc) % p
            else:
                fa = (c * vp[b]) % p; fb = (c * vp[a]) % p
                ga = grad[a]; gb = grad[b]
                if ga and fa:
                    for k, cc in ga.items(): g[k] = (g.get(k, 0) + fa * cc) % p
                if gb and fb:
                    for k, cc in gb.items(): g[k] = (g.get(k, 0) + fb * cc) % p
        else:
            # degree>2 : general product rule
            for idx in range(len(m)):
                fac = c
                for jj, v in enumerate(m):
                    if jj != idx: fac = (fac * vp[v]) % p
                gv = grad[m[idx]]
                if gv and fac:
                    for k, cc in gv.items(): g[k] = (g.get(k, 0) + fac * cc) % p
    return {k: v2 for k, v2 in g.items() if v2}

core = sorted(res)
print("[probe] core root gradient supports:")
for i in core:
    rg = root_grad(i)
    print(f"  eq {i}: residual={res[i]!=0}, |grad|={len(rg)}")
