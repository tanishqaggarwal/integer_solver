#!/usr/bin/env python3
"""Search frame extensions: add atoms to the 7-atom residual frame and recompute max satisfiable rows."""
import sys,os,json,time,itertools
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from full import Full
from fwd import NV
from frame import analyze
F=Full()
d=json.load(open(os.path.join(HERE,'..','best','new_instance_partial_39026.json'))); v=[0]*NV
for k,x in d.items(): v[int(k[2:])]=int(x)
F.init(v)
S0=[j for j in range(F.NA) if F.av[j]]
base=analyze(F,S0)
print('baseline nfail',base['nfail'],flush=True)
E0=set(base['E'])
cand=sorted({j for e in E0 for k,j in F.eqrows[e]} - set(S0))
print('candidates',len(cand),flush=True)
best=base['nfail']
for a in cand:
    t0=time.time()
    r=analyze(F,S0+[a],maxrows=16)
    if r is None:
        print('skip A%d (too many rows)'%a,flush=True); continue
    print('add A%-6d rows=%2d sat=%2d nfail=%2d  t=%.0fs  %s'%(a,len(r['E']),r['nsat'],r['nfail'],time.time()-t0,F.names[a][:60]),flush=True)
    if r['nfail']<best:
        best=r['nfail']; json.dump({'add':a,'frame':r['S'],'E':r['E'],'rows':r['rows'],'sol':{str(k):str(vv) for k,vv in r['sol'].items()}},open(os.path.join(HERE,'frame_better_%d.json'%a),'w'))
        print('  *** BETTER FRAME nfail',best,flush=True)
print('best nfail',best)
