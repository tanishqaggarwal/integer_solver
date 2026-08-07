"""ISD max-satisfy on an ENLARGED 39,026 region (explicit extra atoms)."""
import sys, json, random, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
from gmax import int_solve
P=env.P; Q=(1<<61)-1
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
E=sorted(set(e for a in [22229,22230,35758,35759,35760,35761,35762] for e in L.atom2eq[a]))
A0=set(a for e in E for a in L.eq_atoms[e][2])
EXTRA=[37887,41906,1465,8263,36088,40005,40121,1459,8261,2202,16897,21113,38521,
       39166,40066,40932,29090,2200,21114,32910,1461]
A=A0|set(EXTRA)
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
skipped=len(rows)-len(good)
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Nq=[[x%Q for x in r] for r in N]
d0=[v[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('ENLARGED: atoms=%d knobs=%d eqs=%d skipped=%d nontrivial rows=%d violated=%d %s'%(
      len(A),nk,len(R),skipped,n,len(cur),cur),flush=True)
def infoset(I):
    M=[[Nq[i][j] for j in range(nk)]+[1 if k==t else 0 for t in range(len(I))] for k,i in enumerate(I)]
    r=0
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
        r+=1
    if r!=nk: return None
    Uc=[[M[k][nk+t] for t in range(len(I))] for k in range(nk)]
    return [[sum(Nq[i][j]*Uc[j][t] for j in range(nk))%Q for i in range(n)] for t in range(nk)]
best={}; random.seed(11); t0=time.time(); trials=0
LIM=float(sys.argv[1]) if len(sys.argv)>1 else 600
IDX=list(range(n)); lim=6
while time.time()-t0<LIM:
    trials+=1
    I=random.sample(IDX,nk)
    cws=infoset(I)
    if cws is None: continue
    OUT=[i for i in IDX if i not in I]
    for t in range(nk):
        sup=frozenset(EQ[i] for i in IDX if cws[t][i])
        if 0<len(sup)<=lim: best.setdefault(sup,1)
    for t1 in range(nk):
        for t2 in range(t1+1,nk):
            c1,c2=cws[t1],cws[t2]
            for k in OUT:
                if not c2[k]: continue
                lam=(-c1[k])*pow(c2[k],Q-2,Q)%Q
                sup=frozenset(EQ[i] for i in IDX if (c1[i]+lam*c2[i])%Q)
                if 0<len(sup)<=lim: best.setdefault(sup,2)
print('trials=%d supports<=6 found=%d sizes=%s'%(trials,len(best),
      sorted(collections.Counter(len(s) for s in best).items())),flush=True)
hit=None
for sup in sorted(best,key=len):
    D=set(sup); Z=[i for i in range(n) if EQ[i] not in D]
    x=int_solve([N[i] for i in Z],[B[i] for i in Z],nk)
    if x is not None:
        print('*** INTEGRAL |D|=%d %s'%(len(sup),sorted(sup)),flush=True); hit=(sup,x); break
if hit:
    w=list(v)
    for j,u in enumerate(K): w[u]=hit[1][j]
    av2=L.all_atom_values(w); s2=L.NEQ-len(L.failing_eqs(av2))
    print('SCORE %d'%s2,flush=True)
    if s2>=39026:
        json.dump({str(i):str(w[i]) for i in range(L.NVARS)},
                  open('/home/user/integer_solver/solve_lab/agentA_work/A_enl_%d.json'%s2,'w'))
else:
    print('no integral improvement in the enlarged region',flush=True)
