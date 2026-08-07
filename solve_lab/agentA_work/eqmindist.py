"""Condition (b) done properly at equation level.
(1) RIGOROUS floor: minimum support weight of the code C = colspace(N) equals the minimum
    number of linearly DEPENDENT COLUMNS of a parity check H (basis of the left kernel of
    N).  Exhaustive for small sizes.
(2) Heuristic upper bound: information-set decoding with GREEDY information sets built
    from a shuffled row order (uniform random nk-subsets of rows are almost never
    information sets here -- 0 of 3300 at L=6, which is why eqisd.py returned nothing)."""
import sys, json, time, random, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1]; LEV=int(sys.argv[2]); EXH=int(sys.argv[3]); TL=float(sys.argv[4])
v=L.load(path); fe=L.failing_eqs(L.all_atom_values(v))
A=set(a for e in fe for a in L.eq_atoms[e][2])
for _ in range(LEV):
    R=set()
    for a in A: R|=set(L.atom2eq[a])
    A=set(a for e in R for a in L.eq_atoms[e][2])
K,Rr,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
assert len(good)==len(rows)
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Nq=[[x%Q for x in r] for r in N]
print('lev%d: n=%d rows, nk=%d knobs'%(LEV,n,nk),flush=True)
# ---- parity check H = left kernel of N (mod Q) ----
M=[[Nq[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
r=0
for c in range(nk):
    pr=None
    for i in range(r,n):
        if M[i][c]: pr=i;break
    if pr is None: continue
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
    for i in range(n):
        if i!=r and M[i][c]:
            f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
    r+=1
H=[M[i][nk:] for i in range(r,n)]; m=len(H)
print('rank(N)=%d ; parity check H is %d x %d'%(r,m,n),flush=True)
cols=[[H[t][i] for t in range(m)] for i in range(n)]
zero=[i for i in range(n) if not any(cols[i])]
print('zero columns of H (weight-1 codewords): %d'%len(zero),flush=True)
def rank(vs):
    T=[x[:] for x in vs]; rr=0
    for c in range(m):
        pr=None
        for i in range(rr,len(T)):
            if T[i][c]: pr=i;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],Q-2,Q); T[rr]=[x*inv%Q for x in T[rr]]
        for i in range(len(T)):
            if i!=rr and T[i][c]:
                f=T[i][c]; T[i]=[(a-f*b)%Q for a,b in zip(T[i],T[rr])]
        rr+=1
        if rr==len(T): break
    return rr
t0=time.time(); hit=None
for k in range(2,EXH+1):
    cnt=0
    for S in itertools.combinations(range(n),k):
        cnt+=1
        if rank([cols[i] for i in S])<k: hit=(k,[EQ[i] for i in S]); break
    print('  dependent-column search size %d: %d subsets, found: %s [%.0fs]'%(
        k,cnt,hit[1] if hit else 'NONE',time.time()-t0),flush=True)
    if hit: break
if not hit:
    print('  RIGOROUS: no <= %d dependent columns -> minimum support weight >= %d'%(EXH,EXH+1),flush=True)
# ---- (2) greedy information sets ----
def greedy_infoset(order):
    T=[]; I=[]; piv=[]
    for i in order:
        row=Nq[i][:]
        w=row[:]
        for k2,c in enumerate(piv):
            if w[c]:
                f=w[c]; w=[(a-f*b)%Q for a,b in zip(w,T[k2])]
        nzc=[c for c in range(nk) if w[c]]
        if not nzc: continue
        c0=nzc[0]; inv=pow(w[c0],Q-2,Q); w=[x*inv%Q for x in w]
        for k2 in range(len(T)):
            if T[k2][c0]:
                f=T[k2][c0]; T[k2]=[(a-f*b)%Q for a,b in zip(T[k2],w)]
        T.append(w); piv.append(c0); I.append(i)
        if len(I)==nk: break
    return (I,T,piv) if len(I)==nk else None
random.seed(41); t0=time.time(); tr=0; ok=0; minw=10**9; small=set()
order0=list(range(n))
while time.time()-t0<TL:
    tr+=1
    random.shuffle(order0)
    g=greedy_infoset(order0)
    if g is None: continue
    I,T,piv=g; ok+=1
    Ipos={i:k2 for k2,i in enumerate(I)}
    # codeword with c_I = e_t : u solves N_I u = e_t ; then c = N u
    # T is the RREF of rows I in knob space; build the inverse action
    # easier: solve N_I u = e_t by back-substitution using T/piv on an augmented system
    Aug=[Nq[i][:]+[1 if k2==t else 0 for t in range(nk)] for k2,i in enumerate(I)]
    rr=0; pv=[]
    for c in range(nk):
        pr=None
        for i2 in range(rr,nk):
            if Aug[i2][c]: pr=i2;break
        if pr is None: continue
        Aug[rr],Aug[pr]=Aug[pr],Aug[rr]
        inv=pow(Aug[rr][c],Q-2,Q); Aug[rr]=[x*inv%Q for x in Aug[rr]]
        for i2 in range(nk):
            if i2!=rr and Aug[i2][c]:
                f=Aug[i2][c]; Aug[i2]=[(a-f*b)%Q for a,b in zip(Aug[i2],Aug[rr])]
        pv.append(c); rr+=1
    if rr!=nk: continue
    Uc=[[Aug[k2][nk+t] for t in range(nk)] for k2 in range(nk)]
    for t in range(nk):
        u=[Uc[j][t] for j in range(nk)]
        cw=[sum(Nq[i2][j]*u[j] for j in range(nk))%Q for i2 in range(n)]
        s=frozenset(i2 for i2 in range(n) if cw[i2])
        if s:
            if len(s)<minw:
                minw=len(s)
                print('  greedy-ISD trial %d: support weight %d [%.0fs]'%(tr,minw,time.time()-t0),flush=True)
            if len(s)<=6: small.add(s)
print('greedy-ISD: %d trials, %d usable, lightest support SEEN = %d, supports<=6 found = %d'%(
      tr,ok,minw,len(small)),flush=True)
