import sys, json, math
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, sparse
s=dict(simO.C.base)
v0=E.forward(s); bad0=E.badatoms(v0)
print('bad0',sorted(bad0),'x_15298',v0[15298],'x_7715',v0[7715],'x_34554',v0[34554])
S,cols,nonlin,rounds=simO.closure(v0,bad0,set(),6,8000,verbose=True)
print('knobs',len(S))
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
for a in sorted(bad0):
    reach={f:cols[f][a] for f in S if a in cols[f]}
    lin={f:c for f,c in reach.items() if (f,a) not in nonlin}
    print(f'row {a}: "{H.atoms[a][:80]}" reach={len(reach)} lin={len(lin)} rhs_bits={abs(bad0[a]).bit_length()}')
    for f,c in sorted(lin.items()):
        print(f'   x_{f}: coef({len(str(abs(c)))}d) p|c={c%P==0} c={str(c)[:50]}')
# joint solve on just these two rows
use=sorted(bad0)
rowmap={a:{f:cols[f][a] for f in S if a in cols[f] and (f,a) not in nonlin} for a in use}
sol,msg,_=sparse.solve_sparse([rowmap[a] for a in use],[-bad0[a] for a in use],names=use,verbose=True,maxcore=600,maxcorebits=8_000_000)
print('2-row solve:',msg, 'sol' if sol is not None else 'NONE')
if sol is not None:
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print('EXACT fails',len(ff),'score',39033-len(ff),'bad',sorted(av))
    json.dump({str(k):str(int(x)) for k,x in ns.items()}, open('/home/user/integer_solver/solve_lab/agentO_work/o_2row_seed.json','w'))
