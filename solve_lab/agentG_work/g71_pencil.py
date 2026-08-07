"""The 20 higher-degree equations all lie in the rank-2 pencil spanned by A and B.
Write each as alpha_i*A + beta_i*B and group by direction (alpha:beta).  With (A,B) != 0
the number that can vanish simultaneously is exactly the largest such group."""
import os, sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
W='/home/user/integer_solver/solve_lab/agentG_work/'
D=pickle.load(open(W+'coset_model.pkl','rb')); Lin=pickle.load(open(W+'coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; x0=Lin['x0']
ix={u:i for i,u in enumerate(NB)}
hi=[(i,f) for i,f in non if not isinstance(f,int)]
print('higher-degree equations in the model: %d'%len(hi))
mons=sorted({m for _,f in hi for m in f})
mi={m:t for t,m in enumerate(mons)}
V=[[f.get(m,0) for m in mons] for _,f in hi]
def rank(R):
    R=[r[:] for r in R]; piv=[]; rr=0
    for c in range(len(mons)):
        pr=None
        for t in range(rr,len(R)):
            if R[t][c]%P: pr=t;break
        if pr is None: continue
        R[rr],R[pr]=R[pr],R[rr]
        iv=pow(R[rr][c],-1,P); R[rr]=[x*iv%P for x in R[rr]]
        for t in range(len(R)):
            if t!=rr and R[t][c]%P:
                f=R[t][c]; R[t]=[(x-f*y)%P for x,y in zip(R[t],R[rr])]
        piv.append(c); rr+=1
    return rr,piv,R
rr,piv,RR=rank(V)
print('monomials %d ; SPAN DIMENSION of the higher-degree equations: %d'%(len(mons),rr))
# express each in terms of the first rr independent ones
basis=[]; bidx=[]
cur=[]
for t,(i,f) in enumerate(hi):
    trial=cur+[V[t]]
    if rank(trial)[0]>len(cur): cur=trial; basis.append(V[t]); bidx.append(i)
    if len(cur)==rr: break
print('pencil basis equations: %s'%bidx)
# solve V[t] = sum c_j basis[j]
def solve_coords(v):
    M=[[basis[j][c] for j in range(rr)]+[v[c]] for c in range(len(mons))]
    # gaussian elimination on rr unknowns
    R=[r[:] for r in M]; piv=[]; rr2=0
    for c in range(rr):
        pr=None
        for t in range(rr2,len(R)):
            if R[t][c]%P: pr=t;break
        if pr is None: continue
        R[rr2],R[pr]=R[pr],R[rr2]
        iv=pow(R[rr2][c],-1,P); R[rr2]=[x*iv%P for x in R[rr2]]
        for t in range(len(R)):
            if t!=rr2 and R[t][c]%P:
                fq=R[t][c]; R[t]=[(x-fq*y)%P for x,y in zip(R[t],R[rr2])]
        piv.append(c); rr2+=1
    sol=[0]*rr
    for t,c in enumerate(piv): sol[c]=R[t][rr]%P
    # verify
    for c in range(len(mons)):
        if (sum(basis[j][c]*sol[j] for j in range(rr))-v[c])%P: return None
    return tuple(sol)
grp=collections.defaultdict(list)
for t,(i,f) in enumerate(hi):
    co=solve_coords(V[t])
    if co is None: grp[('NOT-IN-SPAN',)].append(i); continue
    j0=next((j for j in range(rr) if co[j]%P),None)
    if j0 is None: grp[('ZERO',)].append(i); continue
    iv=pow(co[j0],-1,P)
    grp[tuple(c*iv%P for c in co)].append(i)
print('\ndirections in the pencil and the equations carrying them:')
for kdir,v in sorted(grp.items(),key=lambda kv:-len(kv[1])):
    print('   multiplicity %2d : %s'%(len(v),v))
mx=max(len(v) for v in grp.values())
print('\nLARGEST group of proportional higher-degree equations: %d'%mx)
print('=> with the pencil value NONZERO at most %d of the %d can vanish, so at least %d fail.'%(mx,len(hi),len(hi)-mx))
