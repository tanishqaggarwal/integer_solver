#!/usr/bin/env python3
"""Build mod-p Jacobian of the 604 deg-1 constraints w.r.t. support free inputs, then
test consistency of {keep 602 checks=0, fix G1 res, fix G2 res} and solve for delta."""
import sg2_tl as T
import pickle, time, sys
p = T.p
S = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/support.pkl','rb'))
constraints = S['constraints']; rel_order = S['rel_order']
support_free = sorted(S['support_free'])
gate_inputs = S['gate_inputs']
relevant = S['relevant']

vA = T.H.loadd('best_agentA_39022.json')
T.set_base(vA)

# current residuals mod p of the constraints
res = {ai: T.vatom(T.A[ai]['poly']) for ai in constraints}
G1, G2 = 20862, 20864
print("G1 res mod p:", res[G1])
print("G2 res mod p:", res[G2])
homog = [ai for ai in constraints if ai not in (G1,G2)]
print(f"homogeneous checks: {len(homog)}; targets: G1,G2")

# --- precompute descendant-restricted propagation ---
# For efficiency: for a free-input perturbation, propagate through rel_order.
# Build gate ast list restricted (reuse T.gate_ast). Represent tangent of each gate as
# a function using T.val (base) and T.dval.
gate_ast = T.gate_ast
val = T.val; dval = T.dval

# map constraint -> list of (var, is used) for datom; deg-1 so linear: datom = sum coef*dval[var]
# Precompute for each constraint its monomials
cons_terms = {}
for ai in constraints:
    poly = T.A[ai]['poly']
    terms = []  # (coef, tuple-of-vars)
    for m,c in poly.items():
        terms.append((c,m))
    cons_terms[ai] = terms

def datom_fast(ai):
    s = 0
    for c,m in cons_terms[ai]:
        if not m: continue
        # d(prod) = sum_j coef * prod_{k!=j} val[u_k] * dval[u_j]
        L = len(m)
        for j in range(L):
            dv = dval[m[j]]
            if dv == 0: continue
            term = c
            for k in range(L):
                if k==j: term = (term*dv) % p
                else: term = (term*val[m[k]]) % p
            s = (s+term) % p
    return s % p

# Build Jacobian columns. J_rows[i] = dict{free_col: coef}
ci = {ai:i for i,ai in enumerate(constraints)}
J_rows = [dict() for _ in constraints]
t0 = time.time()
for n,f in enumerate(support_free):
    # zero dval on support + relevant
    for v in support_free: dval[v]=0
    for v in relevant: dval[v]=0
    dval[f] = 1
    # propagate through relevant gates in topo order
    for t in rel_order:
        dval[t] = T.devast(gate_ast[t])
    # record datom for each constraint
    for ai in constraints:
        d = datom_fast(ai)
        if d != 0:
            J_rows[ci[ai]][f] = d
    if n % 400 == 0:
        print(f"  col {n}/{len(support_free)}  t={time.time()-t0:.0f}s", flush=True)
print(f"built Jacobian in {time.time()-t0:.0f}s")
# report row nnz stats
nnz = [len(r) for r in J_rows]
print(f"row nnz: min={min(nnz)} max={max(nnz)} mean={sum(nnz)/len(nnz):.1f}")
pickle.dump({'J_rows':J_rows,'constraints':constraints,'support_free':support_free,
             'res':res,'ci':ci}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac.pkl','wb'))
print("saved jac.pkl")
