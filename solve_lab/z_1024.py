import os,sys,json,time,itertools
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
E=FAILS11+RIPPLE16
bits10=[2081,4287,5910,11368,13195,17406,18022,22562,23751,28005]
# free ancestors of E
Eanc=set()
for i in E:
    for v in H.eqvars[i]: Eanc|=H.anc.get(v,{v})
cont=sorted(v for v in Eanc if v in H.freeinp and v not in bits10)
print(f"E={len(E)} eqs, cont free inputs={len(cont)}, bits={len(bits10)}")
# descendant recompute limited to gates feeding E
Egates=set()
for i in E:
    for v in H.eqvars[i]:
        if v not in H.freeinp: Egates.add(v)
# desc within Egates for each input
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    if t in Egates:
        for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid(): return [eval(H.eqcode[i],ns)%p for i in E]
# verify linearity in cont at agentA pattern (2nd difference on a sample)
import random; random.seed(0)
samp=random.sample(cont,min(15,len(cont)))
r0=resid(); ok=True
for v in samp:
    o=H.val[v]; setf(v,o+1); r1=resid(); setf(v,o+2); r2=resid(); setf(v,o)
    for k in range(len(E)):
        if (r2[k]-r1[k]-(r1[k]-r0[k]))%p!=0: ok=False;break
    if not ok:break
print("linear in cont (mod p)?",ok)
# per-pattern feasibility
def close_and_test():
    # close gaps mod p
    setf(7068, H.val[2099]%p); setf(4432, H.val[19964]%p)
    b=resid()
    m=len(E); n=len(cont)
    # jacobian
    cols=[]
    base_cont={v:H.val[v] for v in cont}
    Jt=[[0]*m for _ in cont]  # transpose rows=cont
    for j,v in enumerate(cont):
        o=base_cont[v]; setf(v,o+1); r1=resid(); setf(v,o)
        for k in range(m): Jt[j][k]=(r1[k]-b[k])%p
    # build [J | b] : m x (n+1)
    M=flint.fmpz_mod_mat(m,n+1,ctx)
    for k in range(m):
        for j in range(n): M[k,j]=Jt[j][k]
        M[k,n]=(-b[k])%p
    rJb=M.rank()
    MJ=flint.fmpz_mod_mat(m,n,ctx)
    for k in range(m):
        for j in range(n): MJ[k,j]=Jt[j][k]
    rJ=MJ.rank()
    return rJ,rJb
# test agentA pattern (2081=1 rest0)
for b in bits10: setf(b, 1 if b==2081 else 0)
rJ,rJb=close_and_test()
print(f"agentA pattern (2081=1): rank(J)={rJ}, rank([J|b])={rJb} -> {'FEASIBLE' if rJ==rJb else 'INFEASIBLE (obstruction dim %d)'%(rJb-rJ)}")
