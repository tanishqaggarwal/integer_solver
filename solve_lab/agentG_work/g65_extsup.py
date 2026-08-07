"""Extended departure support: every unknown whose linear footprint lies inside the
region reachable from the deliverable's support (plus the two cheapest unknowns)."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; col=Lin['col']; ix={u:i for i,u in enumerate(NB)}
S0=[642,1329,7068,8731,9118,9413,10903,14623,14853,17325,24548,28730,29854,31339,31864]
U=set()
for u in S0: U|=set(col.get(ix[u],[]))
U |= {56,133,8073}
print('region equations: %d'%len(U))
cand=[]
for c,eqs in col.items():
    if set(eqs)<=U: cand.append((len(eqs),NB[c]))
cand.sort()
extra=[u for _,u in cand if u not in S0]
print('unknowns whose ENTIRE linear footprint lies in the region: %d'%len(cand))
print('  already in the support: %d ; new: %d'%(len(cand)-len(extra),len(extra)))
for w,u in cand:
    print('   x%-6d %d eqs%s'%(u,w,'' if u in S0 else '   <-- NEW'))
SUP=sorted(set(S0)|set(extra)|{22162,30213})
print('\nextended support (%d): %s'%(len(SUP),','.join(map(str,SUP))))
open('/home/user/integer_solver/solve_lab/agentG_work/extsup.txt','w').write(','.join(map(str,SUP)))
