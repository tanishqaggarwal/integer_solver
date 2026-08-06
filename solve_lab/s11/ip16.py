"""IP #16 -- is the RHS p-divisibility reachable using only EXACT-LINEAR variables?
   (ip14 used every variable, including quadratic ones, so its GF(p) 'yes' was not realisable.)"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import atomval, load_raw, deltas
from ip14 import gf_solve
from ip15 import region_of
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); sys.set_int_max_str_digits(400000)
LAB=os.path.join(HERE,'..')
for rel,nm in [(os.path.join(LAB,'best','new_instance_partial_39026.json'),'checkpoint 39026'),
               (os.path.join(HERE,'data','finish3_named.json'),'s11 best 39018'),
               (os.path.join(HERE,'data','closehit2.json'),'closehit2 39005')]:
    v=load_raw(rel)
    AV=[atomval(a,v) for a in range(L.NA)]
    def eqs(e): return sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())
    FAIL=[e for e in range(L.NEQ) if eqs(e)!=0]
    E,cand=region_of(v,AV,FAIL)
    idx={e:i for i,e in enumerate(E)}
    lin=[];allv=[]
    for u in cand:
        d1=deltas(v,AV,u,1)
        if not d1: continue
        allv.append(u)
        d2=deltas(v,AV,u,2)
        if all(d2.get(e,0)==2*d1.get(e,0) for e in set(d1)|set(d2)):
            lin.append((u,d1))
    print(f"=== {nm}: failing={len(FAIL)} region={len(E)} eqs; moving vars={len(allv)}, EXACT-LINEAR={len(lin)}", flush=True)
    if not lin: print("    none"); continue
    M=[[d.get(e,0)%P for u,d in lin] for e in E]
    rhs=[(-eqs(e))%P for e in E]
    t0=time.time(); x=gf_solve(M,rhs,P)
    print(f"    GF(p) solvable with EXACT-LINEAR vars only: {x is not None}  ({time.time()-t0:.0f}s)", flush=True)
