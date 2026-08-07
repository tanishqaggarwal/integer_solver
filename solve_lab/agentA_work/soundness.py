"""SOUNDNESS CHECK of my equation-level formulation.

Question: is my code support vacuous in the way the raw atom-level minimum distance is?
The raw relaxation min_{a != 0} ||M a||_0 <= 1 ranges over ARBITRARY atom vectors a; a
single-equation atom gives weight 1 immediately.

My code is DIFFERENT: C = { N u : u in Q^K } where K are actual VARIABLES and
N[e][j] = d(equation e)/d(knob j).  So a support D means: there is a nonzero KNOB
DIRECTION u with n_e.u = 0 for every e outside D.  Realizability is in the construction.

This script tests the exposure empirically:
 (1) global atom-occupancy distribution and the low-occupancy atoms;
 (2) for every low-occupancy atom IN a window, is D = eqs(a) actually a code support,
     i.e. does rank(N restricted to the rows outside D) drop below #knobs?
 (3) can any knob direction isolate a single atom at all?
"""
import sys, collections, itertools, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
occ=collections.Counter()
for a in range(L.NA): occ[len(L.atom2eq.get(a,{}))]+=1
print('GLOBAL atom occupancy (atoms per #equations):',
      sorted((k,v) for k,v in occ.items() if k<=8),flush=True)
low=[a for a in range(L.NA) if len(L.atom2eq.get(a,{}))<=7]
print('atoms occurring in <= 7 equations: %d ; in exactly 1: %d'%(len(low),occ[1]),flush=True)
print('is a39032 one of them? %s (it occurs in %d equations: %s)'%(
      39032 in low, len(L.atom2eq.get(39032,{})), sorted(L.atom2eq.get(39032,{}))),flush=True)
# does a39032 have a private variable (a variable occurring in no other atom)?
priv=[u for u in L.avars[39032] if len(L.var_atoms[u])==1]
print('a39032 private variables (would make it independently settable): %s'%priv,flush=True)
path='/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
v=L.load(path); fe=L.failing_eqs(L.all_atom_values(v))
for LEV in [int(x) for x in sys.argv[1:]] or [2,6]:
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for _ in range(LEV):
        R=set()
        for a in A: R|=set(L.atom2eq[a])
        A=set(a for e in R for a in L.eq_atoms[e][2])
    K,Rr,rows=build(v,A); nk=len(K)
    NZ=[(e,c,lin) for e,c,lin,hq in rows if not hq and lin]
    n=len(NZ); EQ=[e for e,_,_ in NZ]; E2I={e:i for i,e in enumerate(EQ)}
    N=[[lin.get(j,0)%Q for j in range(nk)] for e,c,lin in NZ]
    def rank(idx):
        M=[N[i][:] for i in idx]; r=0
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
            if r==nk: break
        return r
    full=rank(range(n))
    inA=[a for a in A if len(L.atom2eq.get(a,{}))<=7]
    print('\n--- L=%d : %d atoms, %d rows, %d knobs, rank(N)=%d ---'%(LEV,len(A),n,nk,full),flush=True)
    print('    low-occupancy atoms (<=7 eqs) inside this window: %d'%len(inA),flush=True)
    hits=[]
    for a in inA:
        D=[E2I[e] for e in L.atom2eq[a] if e in E2I]
        if not D: continue
        rr=rank([i for i in range(n) if i not in set(D)])
        if rr<nk: hits.append((a,len(D),sorted(L.atom2eq[a])))
    print('    of those, atoms whose OWN equation set is a code support: %d'%len(hits),flush=True)
    for a,k,es in hits[:6]:
        print('        a%-6d weight %d  eqs %s'%(a,k,es),flush=True)
    # can any knob direction isolate ONE atom?  i.e. is e_a in the image of the knob->atom map?
    # image is a rank-nk lattice inside Z^|A| ; check dimension counting
    print('    knob->atom map: %d knobs -> %d atoms, so the image is a rank-<=%d sublattice'%(
        nk,len(A),nk),flush=True)
    print('    columns of N with the smallest support: %s'%(
        sorted((sum(1 for i in range(n) if N[i][j]),K[j]) for j in range(nk))[:5]),flush=True)
