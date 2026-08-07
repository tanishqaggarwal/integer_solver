"""Light per-boolean scan: which variable does each pin check constrain, and is the
core selector x15298 still live?  Much cheaper than the full reduce."""
import os, sys, time, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G, gGclose
from gsym2 import L, ad, P
SRC='/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v0=L.load(SRC); ad.fwd(v0,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
NB=[u for u in FREE if not gGclose.isbool(u)]
PINS=[1618,688,29539,2423,26731,33796,3576,3578,31670,31672,19297,19299,7930,21617,33929]
def sig(v):
    val,sk=G.build(v,NB,cap=6)
    out={}
    for a in PINS:
        f=G.evalatom(a,val,6)
        if isinstance(f,int): out[a]=('C',f%P)
        else: out[a]=('S',tuple(sorted({NB[k] for m in f for k,_ in m})),G.deg(f))
    return out, v[15298]%P
base,sel0=sig(list(v0))
print('base selector x15298 =',sel0)
for a in PINS: print('  a%-6d %s'%(a,str(base[a])[:100]))
res={}
t0=time.time()
for i,b in enumerate(BOOL):
    w=list(v0); w[b]=1-w[b]; ad.fwd(w,rounds=8)
    try: s,sel=sig(w)
    except Exception as e: res[b]=('ERR',); continue
    diff={a:s[a] for a in PINS if s[a]!=base[a]}
    res[b]=(sel,diff)
    if diff and sel:
        print('bit x%-6d sel=%d CHANGES %s'%(b,sel,{a:(str(v)[:60]) for a,v in list(diff.items())[:6]}),flush=True)
    if i%100==0:
        print('  %d/%d %.0fs'%(i,len(BOOL),time.time()-t0),flush=True)
        pickle.dump({'base':base,'res':res},open('pinscan.pkl','wb'))
pickle.dump({'base':base,'res':res},open('pinscan.pkl','wb'))
nsel=sum(1 for b,r in res.items() if r[0]==sel0 if len(r)>1)
print('bits keeping selector live:',sum(1 for b,r in res.items() if len(r)>1 and r[0]))
print('bits changing any pin:',sum(1 for b,r in res.items() if len(r)>1 and r[1]))
