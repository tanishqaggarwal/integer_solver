"""Price every single-leaf ON configuration through E's simultaneous closure.
   E only ever tried 2 representatives per channel, and always with 1530/1603 ON."""
import sys,json,time,pickle,traceback
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import mcore as M
import mchan as C
import engine as E
BASE=dict(M.load_seed()); BASE[1530]=0; BASE[1603]=0
LEAVES=M.bools()
res={}
best=(10**9,None,None)
t0=time.time()
for i,L in enumerate(LEAVES):
    s=dict(BASE); s[L]=1
    try:
        r=C.simsolve(s)
    except Exception as e:
        res[L]=('ERR',str(e)[:80]); continue
    if r is None: res[L]=('NONE',); continue
    nf,ns,av,v=r
    res[L]=(nf,len(av))
    if nf<best[0]: best=(nf,L,ns)
    if i%16==0:
        print('%3d/%d L=%-6d fails=%-4d score=%-6d best=%d(L=%s) t=%.0fs'%(i,len(LEAVES),L,nf,39033-nf,39033-best[0],best[1],time.time()-t0),flush=True)
pickle.dump(res,open('scan_single.pkl','wb'))
print('BEST fails=%d score=%d leaf=%s'%(best[0],39033-best[0],best[1]))
if best[2] is not None:
    json.dump({str(k):int(x) for k,x in best[2].items()},open('scan_single_best_seed.json','w'))
    v=E.forward(best[2])
    json.dump({("x_%d"%i):v[i] for i in range(E.NV) if v[i]!=0},open('scan_single_best.json','w'))
