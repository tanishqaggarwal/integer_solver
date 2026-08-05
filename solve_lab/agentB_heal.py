#!/usr/bin/env python3
"""Global mod-p heal after the finite control-fix. At the post-fix point, treat the K failing
roots as targets (drive to 0) and all other nonzero rows as constraints (keep 0). Build full
Jacobian, solve the global linear system mod p, report rank/consistency. If consistent, apply and
mod-p Newton-iterate to convergence (handles the quadratic core)."""
import json, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
# finite control fix
x29322 = (env.valp[14853] - env.valp[12186]) % p
x3558 = (env.valp[24908] - env.valp[16742]) % p
env.valp[14853] = (env.valp[14853] - x29322) % p
env.valp[16742] = (env.valp[16742] + x3558) % p
env.forward()

def rref_solve(rows, rhs_map):
    """rows: list of (i, rowdict). rhs_map: i->rhs (target). Solve for delta. Returns
    (delta dict, rank, incons_count, incons_rows)."""
    pivots = {}  # col -> (rowdict, rhs)  normalized
    coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    incons = 0; incons_rows = []
    order = sorted(range(len(rows)), key=lambda k: len(rows[k][1]))
    piv_order = []
    for k in order:
        i, rd0 = rows[k]
        rd = dict(rd0); rhs = rhs_map.get(i, 0)
        while True:
            pcol = None
            for c in rd:
                if c in pivots: pcol = c; break
            if pcol is None: break
            f = rd[pcol]; prow, prhs = pivots[pcol]
            for c, v in prow.items():
                nv = (rd.get(c, 0) - f*v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f*prhs) % p
        if not rd:
            if rhs % p:
                incons += 1; incons_rows.append(i)
            continue
        pc = min(rd, key=lambda c: coldeg.get(c, 0))
        inv = pow(rd[pc], p-2, p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p)
        piv_order.append(pc)
    # back-substitute for a particular solution
    delta = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]
        s = prhs
        for c, v in prow.items():
            if c != pc: s = (s - v * delta.get(c, 0)) % p
        delta[pc] = s % p
    return delta, len(pivots), incons, incons_rows

for it in range(30):
    env.tangent_linear()
    res = env.all_root_residuals()
    if not res:
        print(f"[heal] iter {it}: SOLVED mod p (0 nonzero roots)!"); break
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if g: rows.append((i, g))
    rhs_map = {i: (-res[i]) % p for i in res}
    delta, rank, incons, irows = rref_solve(rows, rhs_map)
    print(f"[heal] iter {it}: failing={len(res)}, rank={rank}, INCONS={incons} {irows[:8]}", flush=True)
    if incons:
        print(f"[heal] iter {it}: linear system INCONSISTENT -> cannot heal globally at this point")
        break
    # apply delta mod p
    for c, d in delta.items():
        env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    res2 = env.all_root_residuals()
    print(f"    after step: failing={len(res2)}")
    if len(res2) >= len(res) and it > 0:
        print("    no progress; stop"); break
