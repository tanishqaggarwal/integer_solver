"""For a given boolean flip set, report the exact inconsistency certificates of the
linear part (which check atoms combine to 0 = nonzero) and the unreachable checks."""
import os, sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym, gGclose
from gsym import *
SRC=os.environ.get('SRC','/home/user/integer_solver/solve_lab/s10/AG_39013.json')
FL=[int(x) for x in sys.argv[1].split(',') if x] if len(sys.argv)>1 and sys.argv[1]!='-' else []
v=L.load(SRC); ad.fwd(v,rounds=6)
for b in FL: v[b]=1-v[b]
ad.fwd(v,rounds=8)
S=gGclose.closure(v); n=len(S)
val=gsym.build(v,S,cap=None,verbose=False)
rows=[];nzc=[]
for a in gsym.check_atoms():
    f=gsym.evalpoly_sym(a,val,n,None)
    if isinstance(f,int):
        if f%P: nzc.append((a,f%P))
    else: rows.append((a,f))
lin=[(a,f) for a,f in rows if gsym.deg(f)==1]
non=[(a,f) for a,f in rows if gsym.deg(f)>1]
m=len(lin)
A=[]
for a,f in lin:
    row=[0]*(n+1+m)
    for mm,c in f.items():
        if sum(mm)==0: row[n]=(-c)%P
        else:
            i=[k for k,e in enumerate(mm) if e][0]; row[i]=c%P
    A.append(row)
for i in range(m): A[i][n+1+i]=1
# rref on first n columns
r=0; piv=[]
for c in range(n):
    pr=None
    for i in range(r,m):
        if A[i][c]%P: pr=i;break
    if pr is None: continue
    A[r],A[pr]=A[pr],A[r]
    iv=pow(A[r][c],-1,P); A[r]=[x*iv%P for x in A[r]]
    for i in range(m):
        if i!=r and A[i][c]%P:
            fq=A[i][c]; A[i]=[(x-fq*y)%P for x,y in zip(A[i],A[r])]
    piv.append(c); r+=1
    if r==m: break
print('flip %s |S|=%d lin=%d rank=%d nonlin=%d nzc=%d'%(FL,n,m,r,len(non),len(nzc)))
nin=0
for i in range(m):
    if all(A[i][c]%P==0 for c in range(n)) and A[i][n]%P:
        nin+=1
        comb=[(lin[j][0],A[i][n+1+j]%P) for j in range(m) if A[i][n+1+j]%P]
        print('  INCONSISTENT cert %d: %d checks -> %s'%(nin,len(comb),
              [(a,('%d'%c if c<10**9 else 'C')) for a,c in comb]))
        print('     eqs per check:',[(a,len(L.atom2eq.get(a,{}))) for a,_ in comb])
if nzc: print('  unreachable nonzero-constant checks:',[(a,len(L.atom2eq.get(a,{}))) for a,_ in nzc])
