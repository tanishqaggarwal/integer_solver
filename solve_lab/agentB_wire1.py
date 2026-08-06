#!/usr/bin/env python3
"""Wire=1 branch analysis. Force all 220 wire members to sign*1 (activating the ~5547 handles),
re-evaluate, and run the EXACT global consistency analysis: active columns, rank(J_sat),
and whether the failing-root residuals are in the column space (J*delta = -residual solvable)."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']

# baseline wire=p
env.set_from_solution(best)
res_p = set(env.all_root_residuals())
print(f"[wire1] wire=p baseline failing: {len(res_p)}")

# wire=1
env.forced = {v: (s % p) for v, s in wire.items()}   # sign*1
env.set_from_solution(best)   # applies forced inside forward()
res1 = env.all_root_residuals()
print(f"[wire1] wire=1 failing: {len(res1)}")
newbreak = sorted(set(res1) - res_p)
healed = sorted(res_p - set(res1))
print(f"[wire1] NEW breaks vs wire=p: {len(newbreak)}: {newbreak}")
print(f"[wire1] previously-failing now healed by wire=1 alone: {len(healed)}: {healed}")

# verify wire members are all 1
bad = [v for v, s in wire.items() if env.valp[v] != (s % p)]
print(f"[wire1] wire members not at sign*1 after forward: {len(bad)}")

# tangent-linear at wire=1: active columns
env.tangent_linear()
active = set()
nnz = 0
sat_rows = []; core_rows = []
for i in range(len(env.root_poly)):
    g = env.root_grad(i)
    if not g: continue
    for c in g: active.add(c)
    nnz += len(g)
    (core_rows if i in res1 else sat_rows).append((i, g))
print(f"[wire1] active columns={len(active)} (was 3036 at wire=p); nnz={nnz}; "
      f"sat rows={len(sat_rows)}, failing rows={len(core_rows)}")

def rref(rows_with_rhs, track=True):
    pivots = {}; coldeg = defaultdict(int)
    for rd, _ in rows_with_rhs:
        for c in rd: coldeg[c] += 1
    incons = 0; ilist = []
    for k in sorted(range(len(rows_with_rhs)), key=lambda k: len(rows_with_rhs[k][0])):
        rd = dict(rows_with_rhs[k][0]); rhs = rows_with_rhs[k][1]
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
            if track: rhs = (rhs - f*prhs) % p
        if not rd:
            if track and rhs % p: incons += 1; ilist.append(k)
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p if track else 0)
    return pivots, len(pivots), incons, ilist

# rank of J_sat (wire=1)
piv_sat, rank_sat, _, _ = rref([(g, 0) for _, g in sat_rows], track=False)
print(f"[wire1] rank(J_sat) = {rank_sat}; nullity over active = {len(active)-rank_sat}")

# full consistency: all failing rows as targets, all sat rows as constraints
res_map = {i: env.root_val(i) for i, _ in core_rows}
rows_all = [(g, 0) for _, g in sat_rows] + [(g, (-res_map[i]) % p) for i, g in core_rows]
piv, rank_all, incons, ilist = rref(rows_all, track=True)
print(f"[wire1] FULL system rank={rank_all}, inconsistent rows={incons}")
print(f"[wire1] => wire=1 linearized system is {'CONSISTENT (SOLVABLE!)' if incons==0 else 'INCONSISTENT'}")

# refined: add failing rows one-by-one to J_sat, count independent vs contradictory
def reduce_against(pivots, rd0, rhs0):
    rd = dict(rd0); rhs = rhs0
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
    return rd, rhs
piv2 = {c: (prow, 0) for c, (prow, _) in piv_sat.items()}
indep=0; incons2=0; absorbed=0
for i, g in core_rows:
    rd, rhs = reduce_against(piv2, g, (-res_map[i]) % p)
    if not rd:
        if rhs % p: incons2 += 1
        else: absorbed += 1
    else:
        pc = min(rd); inv = pow(rd[pc],p-2,p)
        piv2[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); indep += 1
print(f"[wire1] failing rows vs wiring: independent(satisfiable)={indep}, absorbed={absorbed}, "
      f"CONTRADICTORY={incons2} / {len(core_rows)}")
