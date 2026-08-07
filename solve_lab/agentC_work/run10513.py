import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
from plan10513 import plan
W='/home/user/integer_solver/solve_lab/agentC_work/'
BI=json.load(open(W+'bitinfo.json'))
s1=[int(b) for b in BI if BI[b]['side']=='s1']
tb=int(sys.argv[1]); n=int(sys.argv[2]) if len(sys.argv)>2 else 4
best=0
for pa in sorted(s1)[:n]:
    r=plan(tb,pa)
    if r is None: continue
    ctrl,det=r
    t=time.time(); sc,v,nz=closure4(ctrl,detach=det,rounds=30,depth=5)
    print('target x_%d partner x_%d -> %d  nz=%s  %.0fs'%(tb,pa,sc,[(a,len(L.atom2eq.get(a,{}))) for a in nz],time.time()-t),flush=True)
    if sc>best:
        best=sc
        json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open(W+'W_%d_%d_%d.json'%(sc,tb,pa),'w'))
print('BEST',best)
