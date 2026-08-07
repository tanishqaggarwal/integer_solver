#!/usr/bin/env python3
"""Test whether moduli that DIVIDE the handle multipliers M are obstructed too."""
import sys,os,json,re,collections,time,pickle
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine
from modm import solve_modm, ModM
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
pat=re.compile(r'-\((\d+)\*x\d+\)\)$')
Ms=collections.Counter()
for a in E.res:
    m=pat.search(a)
    if m: Ms[int(m.group(1))]+=1
print('distinct handle multipliers M:',len(Ms),'sample',list(Ms)[:8],flush=True)
out=os.path.join(HERE,'modm_results'); os.makedirs(out,exist_ok=True)
tests=[]
for M in list(Ms)[:6]: tests.append((M,'M=%d'%M))
tests.append((7376877,'M_7376877'))
tests.append((7376877*3,'3M'))
for M,tag in tests:
    f=os.path.join(out,'Mmod_%s.json'%tag.replace('=','_'))
    if os.path.exists(f): continue
    t0=time.time()
    v,MM,ok=solve_modm(E,M)
    r=MM.run(v); bad=MM.score(r)
    rec=dict(m=str(M),tag=tag,solved=ok,nonzero_atoms=sum(1 for x in r if x),failing_eqs=len(bad),secs=round(time.time()-t0,1))
    json.dump(rec,open(f,'w'))
    print('modulus %-12s solved=%-5s nz_atoms=%-4d failing_eqs=%-5d t=%.1f'%(str(M)[:12],ok,rec['nonzero_atoms'],len(bad),time.time()-t0),flush=True)
