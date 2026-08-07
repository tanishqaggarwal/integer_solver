"""Symbolize EVERY non-boolean free input (no support closure at all): the maximal
sound model with booleans fixed."""
import os, sys, json, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import gGclose
SRC=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
FL=[int(x) for x in sys.argv[2].split(',') if x] if len(sys.argv)>2 and sys.argv[2]!='-' else []
v=L.load(SRC); ad.fwd(v,rounds=6)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
print('frame %s score %d'%(FL, L.NEQ-len(L.failing_eqs(L.all_atom_values(v)))), flush=True)
FREE=[u for u in range(L.NVARS) if u not in L.definer]
NB=[u for u in FREE if not gGclose.isbool(u)]
print('non-boolean free inputs: %d'%len(NB), flush=True)
t0=time.time()
val,skipped=G.build(v,NB,cap=6,verbose=True)
rows=[];nzc=[]
for a in G.check_atoms():
    f=G.evalatom(a,val,6)
    if isinstance(f,int):
        if f%P: nzc.append((a,f%P))
    else: rows.append((a,f))
print('symbolic checks %d  nonzero-const %d  skipped gates %d  %.1fs'%(len(rows),len(nzc),skipped,time.time()-t0))
print('deg hist',dict(sorted(collections.Counter(G.deg(f) for a,f in rows).items())))
print('nzc:',[(a,len(L.atom2eq.get(a,{}))) for a,_ in nzc][:20])
import pickle; pickle.dump({'syms':NB,'rows':rows,'nzc':nzc},open('allsym_%s.pkl'%('_'.join(map(str,FL)) or 'base'),'wb'))
