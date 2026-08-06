#!/usr/bin/env python3
"""Alternating block Gauss-Seidel over GF(p) for the bilinear wire*handle system, from 39013+wire=1.
Each block step (wire-only, then handle-only) is EXACT linear (other factor fixed) -> no 2nd-order
divergence. Also solves the 3 core handles in closed form each handle step. Goal: all 39033 roots=0
mod p, then Dixon lift."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
handles = set(env.freeinp)
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})

def block_solve(allowed):
    """Solve delta over `allowed` columns to zero as many roots as possible (exact linear in that
    block). Returns delta dict."""
    env.jac_free = wireset  # wire members are columns; handles always columns via freeset
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if not g: continue
        gr = {c: v for c, v in g.items() if c in allowed}
        rows.append((i, gr, (-res[i]) % p if i in res else 0))
    coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []
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
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f*prhs) % p
        if not rd:
            continue  # can't fix with this block; leave for other block
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
    delta = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc:
                dv = delta.get(c,0)
                if dv: s = (s - v*dv) % p
        if s: delta[pc] = s
    return delta

def apply(delta):
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()

print(f"[gs] start (wire=1, 39013): fail={len(env.all_root_residuals())}")
for rnd in range(40):
    t = time.time()
    dH = block_solve(handles); apply(dH); nH = len(env.all_root_residuals())
    dW = block_solve(wireset); apply(dW); nW = len(env.all_root_residuals())
    print(f"  rnd {rnd}: after H fail={nH} (|dH|={len(dH)}), after W fail={nW} (|dW|={len(dW)}), {time.time()-t:.1f}s", flush=True)
    if nW == 0:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced),
                     'freevals': {v: env.valp[v] for v in env.freeset}, 'wireset': list(wireset)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  saved mod-p solution"); break
    if rnd > 3 and nW >= nH:
        print("  no progress"); break
