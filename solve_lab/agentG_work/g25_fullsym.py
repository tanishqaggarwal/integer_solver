"""Symbolize EVERY free input, booleans included (multilinear, b^2=b)."""
import os, sys, json, time, collections, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import gGclose
SRC=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
CAP=int(sys.argv[2]) if len(sys.argv)>2 else 8
v=L.load(SRC); ad.fwd(v,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
print('free %d  boolean %d'%(len(FREE),len(BOOL)),flush=True)
t0=time.time()
try:
    val,skipped=G.build(v,FREE,cap=CAP,boolsyms=BOOL,verbose=True)
except OverflowError as e:
    print('CAP HIT at degree %d'%CAP); sys.exit(1)
rows=[];nzc=[]
for a in G.check_atoms():
    f=G.evalatom(a,val,CAP,boolset=None)
    if isinstance(f,int):
        if f%P: nzc.append((a,f%P))
    else: rows.append((a,f))
print('symbolic checks %d nonzero-const %d skipped %d %.0fs'%(len(rows),len(nzc),skipped,time.time()-t0))
print('deg hist',dict(sorted(collections.Counter(G.deg(f) for a,f in rows).items())))
print('total terms',sum(len(f) for a,f in rows))
pickle.dump({'syms':FREE,'bool':BOOL,'rows':rows,'nzc':nzc},open('fullsym.pkl','wb'))
