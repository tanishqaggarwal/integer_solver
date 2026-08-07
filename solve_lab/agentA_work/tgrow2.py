import sys, json, collections, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from regsolve2 import build, qsolve, pick_knobs
P=env.P
path=sys.argv[1]; STEPS=int(sys.argv[2]) if len(sys.argv)>2 else 25
EQB=int(sys.argv[3]) if len(sys.argv)>3 else 4000
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
print('%s score=%d failing=%d'%(path.split('/')[-1],L.NEQ-len(fe),len(fe)),flush=True)
A=set(a for e in fe for a in L.eq_atoms[e][2])
t0=time.time()
for step in range(STEPS):
    K,R,rows=build(v,A)
    Ks=set(K)
    sol,free,incons,r,piv,mat,eqs=qsolve(rows,len(K))
    nq=sum(1 for x in rows if x[3])
    bad=[(K[j],sol[j]) for j in range(len(K)) if sol[j] is not None and sol[j].denominator!=1]
    print('step%-3d atoms=%-5d knobs=%-4d eqs=%-5d skipped=%-4d rank=%-4d free=%-3d incons=%-3d nonint=%-3d %s [%.0fs]'%(
        step,len(A),len(K),len(R),nq,r,len(free),len(incons),len(bad),
        [(u,'p' if s.denominator==P else (str(s.denominator) if s.denominator<10**9 else 'p*%d'%(s.denominator//P))) for u,s in bad][:8],
        time.time()-t0),flush=True)
    json.dump({'A':sorted(A)},open('/home/user/integer_solver/solve_lab/agentA_work/tg_region.json','w'))
    if incons: print('  INCONSISTENT over Q; stop'); break
    if not bad and not free and nq==0:
        print('*** INTEGRAL SOLUTION OF THE WHOLE REGION ***')
        json.dump({'K':K,'sol':[str(x) for x in sol]},open('/home/user/integer_solver/solve_lab/agentA_work/TG_INTSOL.json','w')); break
    add=set()
    for u,_ in bad:
        for a in L.var_atoms[u]:
            for w in L.avars[a]:
                if w in Ks: continue
                add |= set(L.var_atoms[w])
    add-=A
    if not add: print('  nothing to add; stop'); break
    newe=set(e for a in add for e in L.atom2eq[a])-set(R)
    print('   -> adding %d atoms (+%d eqs)'%(len(add),len(newe)),flush=True)
    if len(R)+len(newe)>EQB: print('   budget exceeded; stop'); break
    A|=add
