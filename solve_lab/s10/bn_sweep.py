"""bn_sweep: measure by construction the cost of turning on each free boolean var.

usage: bn_sweep.py <shard> <nshard> <val>
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

shard=int(sys.argv[1]); nsh=int(sys.argv[2]); VAL=int(sys.argv[3])
bools = B.bools_map()
FREEB = sorted(a for a,(u,c) in bools.items() if u in B.FREESET)
v0 = L.load(B.BEST)
BG = B.broken_gates(v0)
base = B.score(v0)[0]
out=[]
t0=time.time()
for i,a in enumerate(FREEB):
    if i % nsh != shard: continue
    u,c = bools[a]
    if v0[u]==VAL: continue
    w=list(v0); w[u]=VAL
    B.fwdb(w,BG,1)
    s,f,av = B.score(w)
    nz = sum(1 for x in av if x)
    out.append((s,a,u,c,v0[u],nz))
    if len(out)%50==0:
        print(f'  shard{shard} {len(out)} best={max(out)[0]} t={time.time()-t0:.0f}s',flush=True)
out.sort(reverse=True)
json.dump(out, open(os.path.join(HERE,f'bn_sweep_{VAL}_{shard}.json'),'w'))
print(f'shard{shard} done n={len(out)} base={base} TOP:')
for r in out[:15]: print('   score=%d a%d x_%d boolc=%d was=%d nzatoms=%d'%tuple(r))
