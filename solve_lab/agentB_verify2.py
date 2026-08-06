#!/usr/bin/env python3
"""Rigorous independent check of the tangent-linear Jacobian mod p.
For (eq i, free input h): recover P(h)=root_i as polynomial in h via Newton forward differences
at h0..h0+K (exact for deg<=K), compute P'(h0) = sum_{k>=1} (-1)^(k-1) D^k P(h0) / k, compare to
tangent-linear. Uses ONLY forward-eval (independent of the product-rule code)."""
import json, random
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
env.tangent_linear()
res = env.all_root_residuals(); core = sorted(res)
base_free = {v: env.valp[v] for v in env.freeset}

K = 20
invk = [0] + [pow(k, p-2, p) for k in range(1, K+1)]

def P_at(i, h, x):
    env.valp[h] = x % p
    env.forward()
    return env.root_val(i)

def deriv_indep(i, h):
    h0 = base_free[h]
    pts = [P_at(i, h, h0 + t) for t in range(K+1)]
    env.valp[h] = h0; env.forward()  # restore
    # forward differences
    d = pts[:]
    deriv = 0; maxk = 0
    col = pts[:]
    for k in range(1, K+1):
        col = [(col[j+1] - col[j]) % p for j in range(len(col)-1)]
        Dk = col[0]  # D^k P(h0)
        if Dk: maxk = k
        deriv = (deriv + (Dk if k % 2 == 1 else -Dk) * invk[k]) % p
    return deriv, maxk

random.seed(3)
tests = []
for i in core[:8]:
    g = env.root_grad(i)
    for h in list(g)[:4]:
        tests.append((i, h))
# also a few satisfied rows
for i in range(len(env.root_poly)):
    if i not in res:
        g = env.root_grad(i)
        if g and len(g) >= 3:
            for h in list(g)[:2]: tests.append((i, h))
            if len(tests) > 40: break

bad = 0; maxdeg = 0
for i, h in tests:
    g = env.root_grad(i)
    tl = g.get(h, 0) % p
    di, mk = deriv_indep(i, h)
    maxdeg = max(maxdeg, mk)
    if di != tl:
        bad += 1
        print(f"  MISMATCH eq {i} col {h}: indep={di} tl={tl} degree~{mk}")
    env.set_from_solution(best); env.tangent_linear()  # full restore for safety
print(f"[verify2] tested {len(tests)} (eq,col) pairs; mismatches={bad}; max degree seen={maxdeg}")
print("[verify2] tangent-linear Jacobian is", "CORRECT (independently confirmed)" if bad == 0 else "WRONG")
