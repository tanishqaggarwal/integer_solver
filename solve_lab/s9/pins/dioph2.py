"""Proper integer solvability (column-HNF) for the (z,t,h) system."""
import sys, pickle, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

d=pickle.load(open('atoms.pkl','rb')); eq_terms=d['eq_terms']
codes,_=H.load_equations()
FAILS=sorted(H.evaluate(codes,BASE))
D1 = evalpoly(polys[22229],BASE)
D2tot = BASE[4432]-BASE[19964]
rows={}
for e in FAILS:
    m,sq,tl=eq_terms[e]
    co={a:c for c,a in tl}
    if e==8680: rows[e]=(1,0,0,D2tot); continue
    c1=co.get(22229,0); c0=co.get(22230,0); c2=co.get(22231,0)
    rows[e]=(c0-c2, c1*P, -c0*P, -(c1*D1 + c2*D2tot))

def solve_int(M, r):
    """Solve M x = r over Z (M: k x n list of lists). Return x or None."""
    k=len(M); n=len(M[0])
    A=[row[:] for row in M]
    U=[[1 if i==j else 0 for j in range(n)] for i in range(n)]  # column ops accumulate: x = U y
    piv=[]
    row=0; col=0
    while row<k and col<n:
        # reduce columns col..n-1 in this row to a single gcd entry at position col
        while True:
            nz=[j for j in range(col,n) if A[row][j]!=0]
            if len(nz)<=1: break
            nz.sort(key=lambda j: abs(A[row][j]))
            j0=nz[0]
            for j in nz[1:]:
                q=A[row][j]//A[row][j0]
                for i in range(k): A[i][j]-=q*A[i][j0]
                for i in range(n): U[i][j]-=q*U[i][j0]
            # swap j0 into col
            if j0!=col:
                for i in range(k): A[i][col],A[i][j0]=A[i][j0],A[i][col]
                for i in range(n): U[i][col],U[i][j0]=U[i][j0],U[i][col]
        nz=[j for j in range(col,n) if A[row][j]!=0]
        if not nz:
            row+=1; continue
        j0=nz[0]
        if j0!=col:
            for i in range(k): A[i][col],A[i][j0]=A[i][j0],A[i][col]
            for i in range(n): U[i][col],U[i][j0]=U[i][j0],U[i][col]
        piv.append((row,col)); row+=1; col+=1
    # now A is lower-triangular-ish in the pivot positions; forward substitute
    y=[0]*n
    rr=r[:]
    for (i,j) in piv:
        s=rr[i]-sum(A[i][jj]*y[jj] for jj in range(n) if jj!=j)
        if A[i][j]==0 or s%A[i][j]: return None
        y[j]=s//A[i][j]
    for i in range(k):
        if sum(A[i][j]*y[j] for j in range(n))!=rr[i]: return None
    x=[sum(U[i][j]*y[j] for j in range(n)) for i in range(n)]
    for i in range(k):
        if sum(M[i][j]*x[j] for j in range(n))!=r[i]: return None
    return x

if __name__=='__main__':
    results={}
    for k in (4,3,2,1):
        found=[]
        for sub in itertools.combinations(FAILS,k):
            M=[list(rows[e][:3]) for e in sub]; r=[rows[e][3] for e in sub]
            s=solve_int(M,r)
            if s is not None: found.append((sub,s))
        results[k]=found
        print(f'size {k}: {len(found)} solvable subsets')
        for sub,s in found[:12]:
            print('   ',sub,'  z,t,h digit-sizes',[len(str(abs(x))) for x in s])
    pickle.dump(results, open('pins/dioph2.pkl','wb'))
