"""THE decisive test for the canonical basin, mirroring what killed all 38,760 candidates
in the 39,026 region.  Reduce the region system mod p.  Let W = left kernel of N mod p
(dim w), Wb its basis, and g_j = Wb[j].B.  Then the retained rows V\D are mod-p CONSISTENT
iff  g in span{ column_i(Wb) : i in D }.  So the minimum number of equations that must be
violated is exactly the minimum number of columns of Wb needed to represent g."""
import sys, json, time, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P
path=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/mod9118_0.json'
MAXK=int(sys.argv[2]) if len(sys.argv)>2 else 6
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0)%P for j in range(nk)] for e,c,lin in NZ]
B=[(-c)%P for e,c,lin in NZ]
print('%s score=%d : region rows=%d knobs=%d'%(path.split('/')[-1],s0,n,nk),flush=True)
zc=[j for j in range(nk) if all(N[i][j]==0 for i in range(n))]
print('knob columns that vanish mod p (p-quantised handles): %d %s'%(len(zc),[K[j] for j in zc]),flush=True)
# left kernel of N mod p
M=[[N[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
r=0
for c in range(nk):
    pr=None
    for i in range(r,n):
        if M[i][c]: pr=i;break
    if pr is None: continue
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
    for i in range(n):
        if i!=r and M[i][c]:
            f=M[i][c]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r])]
    r+=1
Wb=[M[i][nk:] for i in range(r,n)]
w=len(Wb)
print('rank(N mod p) = %d ; left-kernel dim w = %d'%(r,w),flush=True)
g=[sum(Wb[j][i]*B[i] for i in range(n))%P for j in range(w)]
print('g == 0 (no mod-p obstruction at all) ? %s'%all(x==0 for x in g),flush=True)
if all(x==0 for x in g):
    print('=> the mod-p filter imposes NOTHING here; the basin is decided by integrality only.')
    sys.exit(0)
cols=[[Wb[j][i] for j in range(w)] for i in range(n)]
def solvable(S):
    """is g in span{cols[i] : i in S} ?"""
    T=[cols[i][:] for i in S]+[g[:]]
    ncol=w; rr=0; piv=[]
    for c in range(ncol):
        pr=None
        for i in range(rr,len(S)):
            if T[i][c]: pr=i;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],-1,P); T[rr]=[x*inv%P for x in T[rr]]
        for i in range(len(T)):
            if i!=rr and T[i][c]:
                f=T[i][c]; T[i]=[(a-f*b)%P for a,b in zip(T[i],T[rr])]
        piv.append(c); rr+=1
        if rr==len(S): break
    return all(x==0 for x in T[len(S)])
t0=time.time()
for k in range(1,MAXK+1):
    hit=None; cnt=0
    for S in itertools.combinations(range(n),k):
        cnt+=1
        if solvable(S): hit=[EQ[i] for i in S]; break
    print('  |D|=%d : %d subsets, mod-p-consistent found: %s  [%.0fs]'%(
        k,cnt,hit if hit else 'NONE',time.time()-t0),flush=True)
    if hit:
        json.dump({'D':hit},open('/home/user/integer_solver/solve_lab/agentA_work/modp_hit.json','w'))
        break
else:
    print('CONCLUSION: no drop-set of size <= %d is even mod-p consistent, so every integer'%MAXK)
    print('knob vector in this region violates >= %d equations -> the canonical basin at'%(MAXK+1))
    print('mod9118_0 CANNOT beat 39,026.',flush=True)
