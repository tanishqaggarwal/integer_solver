"""Definitive local rigidity test: over ALL 7,273 free inputs, which move the defects (D1,D2) mod p,
and what do they break?  Perturbations +1 and +p."""
import sys, time, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

freeinp=[x for x in range(NV) if x not in definer]
A0,B0=defects(BASE)
NZ0=set(nz(BASE))
print('free inputs:',len(freeinp),' base defects', A0!=0, B0!=0, ' base residual', sorted(NZ0))

out={}
t0=time.time()
for delta,tag in ((1,'+1'),(P,'+p')):
    rows={}
    for i,f in enumerate(freeinp):
        v=list(BASE)
        try: ripple(v,{f:BASE[f]+delta})
        except Exception: continue
        A,B=defects(v)
        dA,dB=(A-A0)%P,(B-B0)%P
        if dA==0 and dB==0: continue
        broke=sorted(set(nz(v))-NZ0)
        rows[f]=(dA,dB,broke)
    out[tag]=rows
    print(f'\n[{tag}] free inputs that move (D1,D2) mod p: {len(rows)}   [{time.time()-t0:.0f}s]')
    clean=[(f,r) for f,r in rows.items() if not r[2]]
    print(f'   of which break NOTHING new: {len(clean)}')
    for f,r in list(rows.items())[:25]:
        print(f'   x_{f}: dA!=0={r[0]!=0} dB!=0={r[1]!=0} breaks {len(r[2])} atoms {r[2][:6]}')
    for f,r in clean[:20]:
        print(f'   CLEAN MOVER x_{f}: dA={r[0]} dB={r[1]}')
pickle.dump(out, open('pins/sens_defect.pkl','wb'))
