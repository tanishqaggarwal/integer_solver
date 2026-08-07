"""Symbolize only the BOOLEAN free inputs (multilinear, b^2=b), continuous fixed."""
import os, sys, time, collections, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import gGclose
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
CAP=int(sys.argv[1]) if len(sys.argv)>1 else 10
v=L.load(SRC); ad.fwd(v,rounds=6)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
BOOL=[u for u in FREE if gGclose.isbool(u)]
t0=time.time()
try:
    val,skipped=G.build(v,BOOL,cap=CAP,boolsyms=BOOL,verbose=True)
except OverflowError:
    print('CAP %d HIT'%CAP); sys.exit(1)
rows=[];nzc=[]
for a in G.check_atoms():
    f=G.evalatom(a,val,CAP)
    if isinstance(f,int):
        if f%P: nzc.append((a,f%P))
    else: rows.append((a,f))
print('CAP=%d symbolic checks %d nonzero-const %d skipped %d %.0fs'%(CAP,len(rows),len(nzc),skipped,time.time()-t0))
print('deg hist',dict(sorted(collections.Counter(G.deg(f) for a,f in rows).items())))
print('total terms',sum(len(f) for a,f in rows))
pickle.dump({'syms':BOOL,'rows':rows,'nzc':nzc},open('boolsym.pkl','wb'))
