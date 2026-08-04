#!/usr/bin/env python3
"""Iterated mod-p Gauss-Newton. Each step: tangent-linearize; build all nonzero rows with
RHS=-root; RREF (min-coldeg pivot); drop linearly-inconsistent rows THIS step; back-substitute a
particular delta; apply; re-forward. Iterate. Reports convergence of #failing roots.
Modes: 'control' (perturb a var, heal), 'core' (from best), 'corefix' (finite control-fix then Newton)."""
import json, sys, time, random
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
CORE = {2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}

def newton_step(env):
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if g: rows.append((i, g, (-res[i]) % p if i in res else 0))
    coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []; incons = 0
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0, rhs = rows[k]; rd = dict(rd0)
        while True:
            pc = None; best_d = None
            for c in rd:
                if c in pivots:
                    pc = c; break
            if pc is None: break
            f = rd[pc]; prow, prhs = pivots[pc]
            for c, v in prow.items():
                nv = (rd.get(c,0)-f*v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f*prhs) % p
        if not rd:
            if rhs % p: incons += 1
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
    delta = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc:
                dv = delta.get(c, 0)
                if dv: s = (s - v*dv) % p
        if s: delta[pc] = s
    return len(res), len(pivots), incons, delta

def run(env, label, niter=20):
    print(f"=== {label} ===", flush=True)
    for it in range(niter):
        t = time.time()
        nfail, rank, incons, delta = newton_step(env)
        for c, d in delta.items():
            env.valp[c] = (env.valp[c] + d) % p
        env.forward()
        nfail2 = len(env.all_root_residuals())
        ncore = len(env.all_root_residuals().keys() & CORE) if False else sum(1 for i in env.all_root_residuals() if i in CORE)
        print(f"  it {it}: fail {nfail}->{nfail2} (core={ncore}), rank={rank}, incons={incons}, "
              f"|delta|={len(delta)}, {time.time()-t:.1f}s", flush=True)
        if nfail2 == 0:
            print("  *** SOLVED mod p ***"); return True
        if nfail2 >= nfail and it >= 2 and incons > 0:
            print("  plateau"); return False
    return False

mode = sys.argv[1] if len(sys.argv) > 1 else 'control'
if mode == 'control':
    env.set_from_solution(best); env.tangent_linear()
    active = set()
    for i in range(len(env.root_poly)):
        if i in CORE: continue
        for c in env.root_grad(i): active.add(c)
    random.seed(11); h = random.choice(sorted(active))
    env.set_from_solution(best); env.valp[h] = (env.valp[h] + 12345) % p; env.forward()
    print(f"[control] perturbed x_{h}; failing={len(env.all_root_residuals())}")
    run(env, f"iterated-Newton heal of reversible break on x_{h}")
elif mode == 'core':
    env.set_from_solution(best)
    run(env, "iterated-Newton on core from best")
elif mode == 'corefix':
    env.set_from_solution(best)
    x29322 = (env.valp[14853] - env.valp[12186]) % p
    x3558 = (env.valp[24908] - env.valp[16742]) % p
    env.valp[14853] = (env.valp[14853] - x29322) % p
    env.valp[16742] = (env.valp[16742] + x3558) % p
    env.forward()
    print(f"[corefix] after finite fix: failing={len(env.all_root_residuals())}")
    run(env, "iterated-Newton after finite control-fix")
