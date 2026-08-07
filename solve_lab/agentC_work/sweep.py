import sys, json, time, os
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
from plan10513 import plan
W='/home/user/integer_solver/solve_lab/agentC_work/'
BI=json.load(open(W+'bitinfo.json'))
bits=sorted(int(b) for b in BI)
shard=int(sys.argv[1]); nsh=int(sys.argv[2])
part_s1=[b for b in bits if BI[str(b)]['side']=='s1']
part_s2=[b for b in bits if BI[str(b)]['side']=='s2']
res=[]; best=0; t0=time.time()
for i,tb in enumerate(bits):
    if i%nsh!=shard: continue
    partners = part_s2[:2] if BI[str(tb)]['side']=='s1' else part_s1[:2]
    for pa in partners:
        r=plan(tb,pa)
        if r is None: continue
        ctrl,det=r
        try: sc,v,nz=closure4(ctrl,detach=det,rounds=14,depth=4)
        except Exception as e: continue
        res.append((sc,tb,pa,len(nz)))
        if sc>best:
            best=sc; print('BEST %d target x_%d partner x_%d  t=%.0f'%(sc,tb,pa,time.time()-t0),flush=True)
            if sc>=39020:
                json.dump({f'x_{j}':v[j] for j in range(L.NVARS) if v[j]!=0},open(W+'SW_%d_%d_%d.json'%(sc,tb,pa),'w'))
    if len(res)%10==0:
        json.dump(res,open(W+'runs/sweep_%d.json'%shard,'w'))
        print('  %d done best %d %.0fs'%(len(res),best,time.time()-t0),flush=True)
json.dump(res,open(W+'runs/sweep_%d.json'%shard,'w'))
print('DONE',shard,best)
