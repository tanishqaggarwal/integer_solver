"""The 11 failing equations depend on exactly three integer knobs
      A = atom 22229 = x_7068 - x_2099 - 7376877*x_642      (A = D1 + p*t,  t free)
      C = atom 22230 = x_28730 - p*x_9413                   (x_28730 = z free, x_9413 = h free)
      B = atom 22231 = x_4432 - x_19964 - x_28730 = D2tot - z
and on nothing else (x_28730 and x_9413 occur in no other atom; atoms 22229/22230/22231/37887
occur in no equation outside the failing 11).  Solve for the largest subset of the 11 that can be
zeroed simultaneously over Z.
"""
import sys, pickle, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

d=pickle.load(open('atoms.pkl','rb')); eq_terms=d['eq_terms']
codes,_=H.load_equations()
FAILS=sorted(H.evaluate(codes,BASE))
D1 = evalpoly(polys[22229],BASE)
D2tot = BASE[4432]-BASE[19964]
print('D1     =',D1)
print('D2tot  =',D2tot)
print('D1 % p =',D1%P,'\nD2tot%p=',D2tot%P)

# rows: coefficients of (z, t, h) and rhs, for each failing equation
rows={}
for e in FAILS:
    m,sq,tl=eq_terms[e]
    co={a:c for c,a in tl}
    if e==8680:
        rows[e]=(1,0,0, D2tot)            # B = 0  <=>  z = D2tot
        continue
    c1=co.get(22229,0); c0=co.get(22230,0); c2=co.get(22231,0)
    rows[e]=(c0-c2, c1*P, -c0*P, -(c1*D1 + c2*D2tot))
for e in FAILS: print(e, rows[e][:3], 'rhs~%d digits'%len(str(abs(rows[e][3]))))

def solve_sys(eqs):
    """Solve the integer linear system rows[e] . (z,t,h) = rhs  for e in eqs.  Return sol or None."""
    import fractions
    M=[list(rows[e][:3])+[rows[e][3]] for e in eqs]
    n=3
    # integer Gaussian elimination (HNF-style) with column pivoting
    piv=[]; r=0
    cols=list(range(n))
    for c in range(n):
        # find row >= r with nonzero in col c
        k=None
        for i in range(r,len(M)):
            if M[i][c]!=0: k=i; break
        if k is None: continue
        M[r],M[k]=M[k],M[r]
        # eliminate below using integer combinations (gcd based)
        for i in range(r+1,len(M)):
            while M[i][c]!=0:
                q=M[r][c]//M[i][c]
                for j in range(n+1): M[r][j]-=q*M[i][j]
                M[r],M[i]=M[i],M[r]
        piv.append((r,c)); r+=1
        if r==len(M): break
    # consistency of the zero rows
    for i in range(r,len(M)):
        if all(M[i][j]==0 for j in range(n)) and M[i][n]!=0: return None
        if any(M[i][j]!=0 for j in range(n)): return None   # shouldn't happen
    # back substitute: free vars = 0
    sol=[None]*n
    for (rr,cc) in reversed(piv):
        s=M[rr][n]
        for j in range(n):
            if j==cc: continue
            s-= M[rr][j]*(sol[j] if sol[j] is not None else 0)
        if s % M[rr][cc]: return None
        sol[cc]=s//M[rr][cc]
    for j in range(n):
        if sol[j] is None: sol[j]=0
    # verify
    for e in eqs:
        a,b,c,rhs=rows[e]
        if a*sol[0]+b*sol[1]+c*sol[2]!=rhs: return None
    return sol

if __name__=='__main__':
    best=[]
    for k in (3,2,1):
        found=[]
        for sub in itertools.combinations(FAILS,k):
            s=solve_sys(list(sub))
            if s is not None: found.append((sub,s))
        print(f'\nsubsets of size {k} simultaneously zeroable: {len(found)} / {len(list(itertools.combinations(FAILS,k)))}')
        for sub,s in found[:10]:
            print('   ',sub,'-> z,t,h sizes', [len(str(abs(x))) for x in s])
        if found and not best: best=found
    pickle.dump(best, open('pins/dioph_best.pkl','wb'))
