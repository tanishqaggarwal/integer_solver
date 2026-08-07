"""Complete affine model of the 39,026 region: 40 affine rows + the sqrt of a37887."""
import sys, json, collections, math; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from amk_model import build, knobpoly, v0
from agrow import model as amodel
P=env.P
A,K,R=build([37887,41906]); nk=len(K)
kp=knobpoly(37887,K,v0)
VARS=sorted(set(i for m in kp for i in m))
def isqrt_exact(x):
    if x<0: return None
    r=math.isqrt(x)
    return r if r*r==x else None
# q_i from diagonal, sign fixed by cross terms with a reference
q=[0]*(nk+1)   # index nk = constant
diag={}
for i in VARS:
    c=kp.get((i,i),0); s=isqrt_exact(c)
    assert s is not None,(i,c)
    diag[i]=s
c0=kp.get((),0); s0=isqrt_exact(c0); assert s0 is not None
ref=VARS[0]
sign={ref:1}
for i in VARS[1:]:
    cij=kp.get(tuple(sorted((ref,i))),0)
    assert cij % (2*diag[ref]*diag[i])==0 or diag[i]==0, (ref,i)
    sign[i]= 1 if cij>0 else -1
sc = kp.get((ref,),0)
sgn0 = 1 if sc>0 else -1
Q=[0]*nk
for i in VARS: Q[i]=sign[i]*diag[i]
Qc = sgn0*s0
# verify Q^2 == kp
chk=collections.defaultdict(int)
chk[()]=Qc*Qc
for i in VARS:
    chk[(i,)]+=2*Qc*Q[i]
    for j in VARS:
        chk[tuple(sorted((i,j)))]+= Q[i]*Q[j] if i!=j else Q[i]*Q[i]
ok=all(chk.get(m,0)==kp.get(m,0) for m in set(chk)|set(kp))
print('perfect-square check:',ok)
if not ok:
    # try flipping the global sign of the constant
    Qc=-Qc
    chk=collections.defaultdict(int); chk[()]=Qc*Qc
    for i in VARS:
        chk[(i,)]+=2*Qc*Q[i]
        for j in VARS: chk[tuple(sorted((i,j)))]+= Q[i]*Q[j] if i!=j else Q[i]*Q[i]
    ok=all(chk.get(m,0)==kp.get(m,0) for m in set(chk)|set(kp))
    print('after constant sign flip:',ok)
assert ok
d0=[v0[u] for u in K]
print('Q at current knobs =', Qc+sum(Q[i]*d0[i] for i in range(nk)))
print('Q = %s'%(' '.join('%+d*x%d'%(Q[i],K[i]) for i in VARS if Q[i])), ' const digits=%d'%len(str(abs(Qc))))
json.dump({'K':K,'Q':[str(x) for x in Q],'Qc':str(Qc)},open('/home/user/integer_solver/solve_lab/agentA_work/qrow.json','w'))
