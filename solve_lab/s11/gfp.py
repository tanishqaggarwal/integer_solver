"""GF(p) Newton over structural targets."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2, gs
P=L.P
def gauss_solve(A,b,p):
    """min-support solution of A x = b over GF(p); returns list or None"""
    m=len(A); n=len(A[0])
    M=[row[:]+[b[i]] for i,row in enumerate(A)]
    piv=[]; r=0
    for c in range(n):
        pr=None
        for i in range(r,m):
            if M[i][c]%p: pr=i; break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],-1,p)
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%p:
                f=M[i][c]
                M[i]=[(M[i][j]-f*M[r][j])%p for j in range(n+1)]
        piv.append(c); r+=1
        if r==m: break
    for i in range(r,m):
        if M[i][n]%p and all(M[i][j]%p==0 for j in range(n)): return None
    x=[0]*n
    for i,c in enumerate(piv): x[c]=M[i][n]%p
    return x

TARGETS=[('x3719',   lambda v: v[3719]),
         ('x25118',  lambda v: v[25118]),
         ('x25614',  lambda v: v[25614]),
         ('x34220',  lambda v: v[34220]),
         ('n-gap',   lambda v: v[12186]-v[1308]),
         ('m-gap',   lambda v: v[24908]-v[19083])]
CTRL=[14515,19750,5096,21589,33708,31339,29261,26489,8060,19450,3473,8971,5616,245]
def resid(th):
    v=gs.state(th)
    return v,[f(v)%P for _,f in TARGETS]
if __name__=='__main__':
    theta={c:0 for c in CTRL}
    t0=time.time()
    for it in range(12):
        v,r=resid(theta)
        print(f"it{it}: residuals nonzero = {[TARGETS[i][0] for i,x in enumerate(r) if x]}  ({time.time()-t0:.0f}s)")
        if not any(r): break
        J=[[0]*len(CTRL) for _ in TARGETS]
        for j,c in enumerate(CTRL):
            th=dict(theta); th[c]=theta[c]+1
            _,r1=resid(th)
            for i in range(len(TARGETS)): J[i][j]=(r1[i]-r[i])%P
        d=gauss_solve(J,[(-x)%P for x in r],P)
        if d is None: print("  GF(p) system INCONSISTENT"); break
        for j,c in enumerate(CTRL): theta[c]=(theta[c]+d[j])%P
    v,r=resid(theta)
    print("all targets zero:", not any(r))
    json.dump({str(k):val for k,val in theta.items()}, open('theta_gfp.json','w'))
