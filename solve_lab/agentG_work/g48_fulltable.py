"""Enumerate the point TABLE: for every boolean free input, which of P1,P2,P3 it moves
and to which secp256k1 point.  Then search for i,j,k with P1(i)+P2(j) = P3(k)."""
import os, sys, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g46_table as T, gpt, gGclose
from gsym2 import L, ad, P
OUT='/home/user/integer_solver/solve_lab/agentG_work/fulltable.pkl'
LO=int(sys.argv[1]) if len(sys.argv)>1 else 0
HI=int(sys.argv[2]) if len(sys.argv)>2 else 100000
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
base=T.frame([])
B=base['pts']
print('bits=%d  base P1,P2,P3 on curve=%s  D=%s'%(len(BOOL),base['oncurve'],base['D']),flush=True)
res={'base':base}
t0=time.time()
for i,b in enumerate(BOOL[LO:HI]):
    try: o=T.frame([b])
    except Exception as e:
        res[b]=('ERR',str(e)[:40]); continue
    rec={'nres':o.get('nres'),'ninc':o.get('ninc'),'nzc':o.get('nzc'),
         'core':o.get('core'),'oncurve':o.get('oncurve'),'pts':o.get('pts'),
         'label':o.get('label'),'D':o.get('D')}
    res[b]=rec
    if o.get('pts'):
        moved=[k for k in range(3) if o['pts'][k]!=B[k]]
        rec['moved']=moved
        if o.get('D') is None and o.get('oncurve') and not all(o['oncurve']):
            pass
        if o.get('D') is not None and o['D'] is None:
            print('*** D = O (P3 = P1+P2) at bit x%d ***'%b,flush=True)
    if i%20==0:
        print('  %d/%d %.0fs  x%-6d moved=%s oncurve=%s nres=%s'%(LO+i,min(HI,len(BOOL)),time.time()-t0,b,rec.get('moved'),rec.get('oncurve'),rec.get('nres')),flush=True)
        pickle.dump(res,open(OUT,'wb'))
pickle.dump(res,open(OUT,'wb'))
print('done %.0fs'%(time.time()-t0))
