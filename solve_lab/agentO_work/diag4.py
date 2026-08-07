import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
s=dict(simO.C.base); s[22492]=1
v0=E.forward(s)
def show(u):
    d=H.definer[u]
    print(f'x_{u}: def={"FREE" if d is None else H.atoms[d[0]][:120]} val_bits={abs(v0[u]).bit_length()} val0={v0[u]==0}')
for u in [15298,11150,25739,29804,37758,35605,22162,30213,13682,18956,32237,6083,14853,36433,36990,19239,24530,5647,24908]:
    show(u)
print('--- atoms')
for a in [20649,20652,32148,28647,20215,20212]:
    print(a, H.atoms[a][:150])
