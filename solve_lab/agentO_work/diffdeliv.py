import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, full11 as F, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
FR=[u for u in range(E.NV) if E.definer[u] is None]
nzd={u:vd[u] for u in FR if vd[u]!=0}
print('deliverable nonzero free vars:',len(nzd))
r=F.solve_pair(24601,2081,verbose=False)
n,ns,av=r
v=E.forward(ns)
nzp={u:v[u] for u in FR if v[u]!=0}
print('full11 pair nonzero free vars:',len(nzp))
print('--- in deliverable only:')
for u in sorted(set(nzd)-set(nzp)): print(f'  x_{u} = {str(nzd[u])[:50]} ({abs(nzd[u]).bit_length()}b)')
print('--- in full11 only:')
for u in sorted(set(nzp)-set(nzd)): print(f'  x_{u} = {str(nzp[u])[:50]} ({abs(nzp[u]).bit_length()}b)')
print('--- both but differ:')
for u in sorted(set(nzp)&set(nzd)):
    if nzp[u]!=nzd[u]: print(f'  x_{u}: deliv {str(nzd[u])[:40]}({abs(nzd[u]).bit_length()}b)  pair {str(nzp[u])[:40]}({abs(nzp[u]).bit_length()}b)')
print('--- same:', sum(1 for u in set(nzp)&set(nzd) if nzp[u]==nzd[u]))
