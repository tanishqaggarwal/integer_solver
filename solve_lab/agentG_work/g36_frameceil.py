import os, sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gframe
from gsym2 import L, ad, P
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
v0=L.load(SRC); ad.fwd(v0,rounds=6)
for arg in sys.argv[1:]:
    parts=arg.split('/')
    FL=[int(x) for x in parts[0].split(',') if x] if parts[0] not in ('-','') else []
    DET=[int(x) for x in parts[1].split(',') if x] if len(parts)>1 and parts[1] not in ('-','') else []
    v=list(v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    t0=time.time()
    nfail,inc,res=gframe.ceiling(v,DET)
    print('frame flips=%s detach=%s -> exact mod-p ceiling %d (fail %d: inc %d res %d)  %.0fs'
          %(FL,DET,L.NEQ-nfail,nfail,len(inc),len(res),time.time()-t0),flush=True)
