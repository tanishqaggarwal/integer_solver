#!/usr/bin/env python3
"""Build the full sparse mod-p root-Jacobian at the best solution and measure structure."""
import json, time, sys
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

t0 = time.time()
data = load()
env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
env.set_from_solution(best)
env.tangent_linear()
res = env.all_root_residuals()
neq = len(env.root_poly)
print(f"[jac] setup+forward+tangent in {time.time()-t0:.1f}s; nonzero roots={len(res)}")

t1 = time.time()
rows = []          # list of (i, dict col->val)
nnz_total = 0
coldeg = defaultdict(int)
nz_rows = 0
for i in range(neq):
    g = env.root_grad(i)
    if g:
        rows.append((i, g))
        nz_rows += 1
        nnz_total += len(g)
        for c in g: coldeg[c] += 1
print(f"[jac] built {nz_rows} nonzero rows / {neq} eqs in {time.time()-t1:.1f}s")
print(f"[jac] total nonzeros={nnz_total}, mean nnz/nonzero-row={nnz_total/max(nz_rows,1):.1f}")
print(f"[jac] distinct free-input columns used={len(coldeg)} / {len(env.freeinp)}")
# row nnz distribution
import statistics
nnzs = [len(g) for _, g in rows]
nnzs.sort()
print(f"[jac] row nnz: min={nnzs[0]} med={nnzs[len(nnzs)//2]} p90={nnzs[int(len(nnzs)*0.9)]} max={nnzs[-1]}")
cds = sorted(coldeg.values())
print(f"[jac] col deg: min={cds[0]} med={cds[len(cds)//2]} p90={cds[int(len(cds)*0.9)]} max={cds[-1]}")
# how many core rows and their columns
core = sorted(res)
print(f"[jac] core rows (nonzero residual): {len(core)}")
# save rows for downstream solvers
import pickle
SC = '/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_jac.pkl'
with open(SC, 'wb') as f:
    pickle.dump({'rows': rows, 'res': res, 'freeinp': env.freeinp}, f)
print(f"[jac] saved {SC} in total {time.time()-t0:.1f}s")
