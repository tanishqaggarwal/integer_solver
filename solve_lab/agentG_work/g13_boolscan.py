import os, sys, json, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gred
from gsym import *
src='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
SYMS=json.load(open('closed_nonbool.json'))
v0=L.load(src); ad.fwd(v0,rounds=6)
# all boolean free inputs
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[]
for u in FREE:
    for a in L.var_atoms[u]:
        pl=L.polys[a]
        if len(pl)==2 and (u,) in pl and (u,u) in pl and pl[(u,)]==-pl[(u,u)]:
            BOOL.append(u); break
print('boolean free inputs:', len(BOOL), flush=True)
base=gred.reduce_state(list(v0),SYMS)
bs=gred.residual_signature(base)
print('base: rank %d free %d nzc %d residual %s'%(base['rank'],base['nfree'],len(base['nzc']),[(a,str(g)[:12]) for a,g in bs]), flush=True)
out={}
t0=time.time()
for i,b in enumerate(BOOL):
    w=list(v0); w[b]=1-w[b]
    ad.fwd(w,rounds=8)
    try:
        r=gred.reduce_state(w,SYMS)
    except Exception as e:
        out[b]=('ERR',str(e)[:40]); continue
    sig=gred.residual_signature(r)
    out[b]=dict(rank=r['rank'],nfree=r['nfree'],nzc=len(r['nzc']),ninc=r['ninc'],
                nres=len([1 for a,g in sig]), sig=[(a,(g if isinstance(g,int) else g)) for a,g in sig])
    if len(sig)==0 and len(r['nzc'])==0:
        print('*** ZERO RESIDUAL at bit x%d'%b, flush=True)
    if i%100==0: print('  %d/%d  %.0fs'%(i,len(BOOL),time.time()-t0), flush=True)
pickle.dump({'base':bs,'out':out}, open('boolscan.pkl','wb'))
# summarize
import collections
cnt=collections.Counter()
for b,r in out.items():
    if isinstance(r,tuple): cnt['ERR']+=1; continue
    key=(r['nres'],r['nzc'],r['ninc'])
    cnt[key]+=1
print('summary (nres,nzc,ninc):',cnt.most_common(20))
best=[(b,r) for b,r in out.items() if not isinstance(r,tuple) and r['nres']<len(bs)]
print('bits reducing residual count:',[(b,r['nres']) for b,r in best][:40])
