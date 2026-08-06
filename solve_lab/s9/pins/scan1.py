"""Per-bit effect: turn bit b on (with its pins closed) from the best partial; record residual+defects."""
import sys, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

BASEA, BASEB = defects(BASE)
res = {}
t0=time.time()
for i,b in enumerate(BITS):
    if BASE[b]: continue
    v=list(BASE)
    s={b:1}
    for pn in bitpins[b]:
        s[pn['B']]=pn['HUGE']; s[pn['h']]=0
    ripple(v,s)
    n=nz(v); A,B=defects(v)
    res[b]=dict(nz=n, dA=(A-BASEA)%P, dB=(B-BASEB)%P, A=A, B=B)
    if i%40==0: print(i,b,len(n),f'{time.time()-t0:.0f}s',flush=True)
pickle.dump(res, open('pins/scan1.pkl','wb'))
import collections
print('done',f'{time.time()-t0:.0f}s')
print('residual-size histogram:', collections.Counter(len(r['nz']) for r in res.values()))
small=sorted(res.items(), key=lambda kv: len(kv[1]['nz']))[:20]
for b,r in small: print(b, len(r['nz']), r['nz'][:12], 'dA!=0',r['dA']!=0,'dB!=0',r['dB']!=0)
