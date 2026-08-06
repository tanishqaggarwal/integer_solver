#!/usr/bin/env python3
"""Build mod-p Jacobian of ALL 4191 sensitive atoms w.r.t support free inputs (optimized),
then solve {keep all sensitive checks=0 (first order), fix G1,G2 res}."""
import sg2_tl as T
import pickle, time
from collections import defaultdict
p = T.p
S = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/support2.pkl','rb'))
constraints = S['constraints']; rel_order = S['rel_order']
support_free = sorted(S['support_free']); relevant = S['relevant']
idx = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/idx.pkl','rb'))
var_atoms = idx['var_atoms']
G1,G2 = 20862,20864
consset = set(constraints)
# var -> constraints containing it
var_cons = defaultdict(list)
for ai in constraints:
    for v in T.L.atom_vars(T.A[ai]['poly']):
        var_cons[v].append(ai)

vA = T.H.loadd('best_agentA_39022.json'); T.set_base(vA)
res = {ai: T.vatom(T.A[ai]['poly']) for ai in constraints}
val=T.val; dval=T.dval; gate_ast=T.gate_ast

# precompute constraint terms
cons_terms = {ai:[(c,m) for m,c in T.A[ai]['poly'].items()] for ai in constraints}
def datom_fast(ai):
    s=0
    for c,m in cons_terms[ai]:
        Lm=len(m)
        if Lm==0: continue
        for j in range(Lm):
            dv=dval[m[j]]
            if dv==0: continue
            term=c
            for k in range(Lm):
                term = (term*(dv if k==j else val[m[k]]))%p
            s=(s+term)%p
    return s%p

ci={ai:i for i,ai in enumerate(constraints)}
J_rows=[dict() for _ in constraints]
t0=time.time()
allsup=support_free
relist=rel_order
for n,f in enumerate(allsup):
    dval[f]=1
    for t in relist:
        dval[t]=T.devast(gate_ast[t])
    # hot vars
    hot=[f]+[t for t in relist if dval[t]]
    cand=set()
    for v in hot: cand.update(var_cons.get(v,()))
    for ai in cand:
        d=datom_fast(ai)
        if d: J_rows[ci[ai]][f]=d
    # reset
    dval[f]=0
    for t in relist: dval[t]=0
    if n%800==0: print(f"  col {n}/{len(allsup)} t={time.time()-t0:.0f}s",flush=True)
print(f"Jacobian built {time.time()-t0:.0f}s")
nnz=[len(r) for r in J_rows]
print(f"row nnz min={min(nnz)} max={max(nnz)} mean={sum(nnz)/len(nnz):.1f}; rows w/ 0 nnz={sum(1 for x in nnz if x==0)}")
pickle.dump({'J_rows':J_rows,'constraints':constraints,'support_free':support_free,'res':res,'ci':ci},
            open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac2.pkl','wb'))
print("saved jac2.pkl")
