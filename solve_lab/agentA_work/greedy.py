"""Greedy enlargement of the atom set A*: at each step add the external-atom set of
the region variable that adds the fewest new equations.  Report knobs/rank/denominators."""
import sys, json, collections, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from amk_model import build, knobpoly, v0, av0, A0
from densolve import qsolve
from agrow import model
P=env.P
BUDGET=int(sys.argv[1]) if len(sys.argv)>1 else 400
A=set(A0)
hist=[]
t0=time.time()
for step in range(200):
    A2,K,R,rows,QUAD=model(sorted(A-set(A0)))
    sol,free,incons,r=qsolve(rows,len(K))
    nq=sum(1 for x in rows if x[3])
    dens=collections.Counter()
    bad=[]
    if not incons:
        for j,u in enumerate(K):
            if sol[j] is None: continue
            d=sol[j].denominator
            if d==1: continue
            bad.append((u,'p' if d==P else ('p*%d'%(d//P) if d%P==0 else str(d))))
    print('step%-3d atoms=%-5d knobs=%-4d eqs=%-5d quadrows=%-3d rank=%-4d free=%-3d incons=%-3d nonint=%-3d %s  [%.0fs]'%(
        step,len(A2),len(K),len(R),nq,r,len(free),len(incons),len(bad),bad[:8],time.time()-t0),flush=True)
    if not incons and not bad and not free:
        print('*** FULLY INTEGRAL SOLUTION FOUND ***'); json.dump(
            {'K':K,'sol':[str(x) for x in sol]},open('/home/user/integer_solver/solve_lab/agentA_work/INTSOL.json','w')); break
    if len(R)>BUDGET: print('budget reached'); break
    # choose the next variable to free: the one whose external atoms add fewest equations
    Aset=set(A2); best=None
    Vs=set(u for a in A2 for u in L.avars[a])
    for u in Vs:
        ext=[a for a in L.var_atoms[u] if a not in Aset]
        if not ext: continue
        newe=set()
        for a in ext: newe |= set(L.atom2eq[a])-set(R)
        cand=(len(newe),len(ext),u,ext)
        if best is None or cand[:3]<best[:3]: best=cand
    if best is None: print('closed'); break
    print('   -> freeing x%d by adding atoms %s (+%d eqs)'%(best[2],best[3][:6],best[0]),flush=True)
    A |= set(best[3])
