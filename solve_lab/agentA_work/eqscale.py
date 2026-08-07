"""How far can the FULL pipeline be pushed?  Time the mod-p reduction per level."""
import sys, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fe=L.failing_eqs(L.all_atom_values(v))
A=set(a for e in fe for a in L.eq_atoms[e][2])
for lev in range(0,int(sys.argv[1])+1):
    if lev:
        R=set()
        for a in A: R|=set(L.atom2eq[a])
        A=set(a for e in R for a in L.eq_atoms[e][2])
    if lev not in (6,8,10,12,14,16): continue
    t0=time.time()
    K,Rr,rows=build(v,A); nk=len(K)
    good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    NZ=[(e,c,lin) for e,c,lin in good if lin]
    n=len(NZ)
    N=[[lin.get(j,0)%P for j in range(nk)] for e,c,lin in NZ]
    tb=time.time()-t0
    t0=time.time()
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
    print('L%-3d n=%-5d nk=%-5d rank_modp=%-5d w=%-5d build=%.1fs kernel=%.1fs'%(
        lev,n,nk,r,n-r,tb,time.time()-t0),flush=True)
