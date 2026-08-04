#!/usr/bin/env python3
"""Verify the tangent-linear Jacobian mod p against finite differences (symmetric, exact for
deg<=2). Also refined rank decomposition: rank(J_sat), and how many of the 20 core conditions
are linearly consistent with the wiring."""
import json, time, random
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
env.tangent_linear()
res = env.all_root_residuals()
core = sorted(res)

# ---- finite-difference verification (mod p) ----
def root_val_mod(i):
    return env.root_val(i)
def perturb_free(h, d):
    env.valp[h] = (env.valp[h] + d) % p
    env.forward()
random.seed(0)
bad = 0; checked = 0
for i in core[:6] + [r[0] for r in [(rr,) for rr in []]]:
    g = env.root_grad(i)
    cols = list(g)[:5]
    for h in cols:
        base = env.valp[h]
        perturb_free(h, 1); rp = root_val_mod(i)
        perturb_free(h, -2); rm = root_val_mod(i)
        perturb_free(h, 1)  # restore
        deriv_fd = ((rp - rm) * pow(2, p-2, p)) % p
        if deriv_fd != g[h] % p:
            bad += 1
            print(f"  MISMATCH eq {i} col {h}: fd={deriv_fd} tl={g[h]%p}")
        checked += 1
print(f"[verify] finite-diff check: {checked} (col,eq) pairs, mismatches={bad}")

# ---- refined rank: build rows, separate satisfied vs core ----
env.set_from_solution(best); env.tangent_linear()
sat_rows = []; core_rows = []
for i in range(len(env.root_poly)):
    g = env.root_grad(i)
    if not g: continue
    if i in res: core_rows.append((i, g, (-res[i]) % p))
    else: sat_rows.append((i, g, 0))
print(f"[verify] satisfied nonzero rows={len(sat_rows)}, core rows={len(core_rows)}")

def rref(rows_with_rhs, track_rhs=True):
    """rows_with_rhs: list of (rowdict, rhs). Returns (pivots dict col->(prow,prhs), rank, incons_count)."""
    pivots = {}; coldeg = defaultdict(int)
    for rd, _ in rows_with_rhs:
        for c in rd: coldeg[c] += 1
    incons = 0; npiv = 0
    order = sorted(range(len(rows_with_rhs)), key=lambda k: len(rows_with_rhs[k][0]))
    for k in order:
        rd, rhs = rows_with_rhs[k]; rd = dict(rd)
        while True:
            pcol = None
            for c in rd:
                if c in pivots: pcol = c; break
            if pcol is None: break
            f = rd[pcol]; prow, prhs = pivots[pcol]
            for c, v in prow.items():
                nv = (rd.get(c, 0) - f * v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            if track_rhs: rhs = (rhs - f * prhs) % p
        if not rd:
            if track_rhs and rhs % p: incons += 1
            continue
        pc = min(rd, key=lambda c: coldeg.get(c, 0))
        inv = pow(rd[pc], p-2, p)
        pivots[pc] = ({c: (v*inv) % p for c, v in rd.items()}, (rhs*inv) % p if track_rhs else 0)
        npiv += 1
    return pivots, npiv, incons

# rank of satisfied constraints alone
piv_sat, rank_sat, _ = rref([(g, 0) for _, g, _ in sat_rows], track_rhs=False)
print(f"[verify] rank(J_satisfied) = {rank_sat}  (cols used among satisfied)")

# now reduce each core row's LHS against J_sat pivots: does it add a new pivot (independent) or not?
def reduce_against(pivots, rd0, rhs0):
    rd = dict(rd0); rhs = rhs0
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
    return rd, rhs
# reduce core rows one by one, accumulating pivots (so we detect dependent core rows)
piv = dict(piv_sat)
core_indep = 0; core_incons = 0; core_absorbed = 0
for i, g, rhs in core_rows:
    rd, r = reduce_against(piv, g, rhs)
    if not rd:
        if r % p: core_incons += 1
        else: core_absorbed += 1
    else:
        pc = min(rd, key=lambda c: len(str(c)))
        inv = pow(rd[pc], p-2, p)
        piv[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (r*inv)%p)
        core_indep += 1
print(f"[verify] core rows: independent(add new pivot, satisfiable)={core_indep}, "
      f"absorbed(0=0)={core_absorbed}, INCONSISTENT(0=nonzero)={core_incons}")
print(f"[verify] => {core_incons} of 20 core conditions are LINEARLY UNREACHABLE from best (mod p)")
