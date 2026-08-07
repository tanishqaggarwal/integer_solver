#!/usr/bin/env python3
"""AUDIT T6b: structural comparison of A's 9-knob L=0 window against the 24-knob one.
rank over Q, Q-consistency, mod-p rank/consistency, and the lightest mod-p-admissible
violated set (agent A's condition (a)) -- computed with the exclusion removed."""
import os,sys,pickle,json,collections,itertools,time
from fractions import Fraction
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
from circ2 import vars_of
from fwd import compile_node
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
avars=[frozenset(vars_of(atoms[a])) for a in names]
v2a=collections.defaultdict(set)
for i,vs in enumerate(avars):
    for u in vs: v2a[u].add(i)
eq_atoms=[frozenset(idx[a] for k,a in row) for row in eqrows]
a2e=collections.defaultdict(set)
for e,s in enumerate(eq_atoms):
    for a in s: a2e[a].add(e)
NV=38748
asg=json.load(open(os.path.join(T,'..','best','new_instance_partial_39026.json')))
V0=[0]*NV
for k,val in asg.items(): V0[int(k[2:])]=int(val)
LEV=int(sys.argv[1]) if len(sys.argv)>1 else 0
FAIL=[12231,12270,12350,14584,18673,22044,29125]
A=set()
for e in FAIL: A|=eq_atoms[e]
for L in range(LEV):
    Rr=set()
    for a in A: Rr|=a2e[a]
    A2=set(A)
    for e in Rr: A2|=eq_atoms[e]
    A=A2
R=sorted({e for a in A for e in a2e[a]})
Vw=set()
for a in A: Vw|=avars[a]
K=sorted(u for u in Vw if v2a[u]<=A)
need=sorted({a for e in R for a in eq_atoms[e]})
sub=[names[i] for i in need]; pos={a:i for i,a in enumerate(need)}
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in sub)+']','<at>','exec')
rr=[0]*len(need)
def cores(v):
    exec(prog,{'v':v,'r':rr,'__builtins__':{}})
    return [sum(k*rr[pos[idx[a]]] for k,a in eqrows[e]) for e in R]
c0=cores(V0)
N=[[0]*len(K) for _ in R]
for j,u in enumerate(K):
    v=list(V0); v[u]=V0[u]+1; c1=cores(v)
    for i in range(len(R)): N[i][j]=c1[i]-c0[i]
B=[-x for x in c0]
def rank_mod(M,n,q):
    M=[[x%q for x in r] for r in M]; piv=0
    for c in range(n):
        k=None
        for i in range(piv,len(M)):
            if M[i][c]%q: k=i;break
        if k is None: continue
        M[piv],M[k]=M[k],M[piv]; inv=pow(M[piv][c],q-2,q)
        M[piv]=[(x*inv)%q for x in M[piv]]
        for i in range(len(M)):
            if i!=piv and M[i][c]%q:
                f=M[i][c]; M[i]=[(M[i][j]-f*M[piv][j])%q for j in range(n)]
        piv+=1
    return piv
def rank_q(M,n):
    M=[[Fraction(x) for x in r] for r in M]; piv=0
    for c in range(n):
        k=None
        for i in range(piv,len(M)):
            if M[i][c]!=0: k=i;break
        if k is None: continue
        M[piv],M[k]=M[k],M[piv]; pv=M[piv][c]
        M[piv]=[x/pv for x in M[piv]]
        for i in range(len(M)):
            if i!=piv and M[i][c]!=0:
                f=M[i][c]; M[i]=[M[i][j]-f*M[piv][j] for j in range(n)]
        piv+=1
    return piv
n=len(K)
Naug=[N[i]+[B[i]] for i in range(len(R))]
rQ=rank_q(N,n); rQa=rank_q(Naug,n+1)
rP=rank_mod(N,n,p); rPa=rank_mod(Naug,n+1,p)
print('L=%d  rows=%d knobs=%d'%(LEV,len(R),n))
print('  rank_Q(N)=%d  rank_Q([N|B])=%d  -> Q-consistent: %s'%(rQ,rQa,rQ==rQa))
print('  rank_p(N)=%d  rank_p([N|B])=%d  -> mod-p consistent: %s'%(rP,rPa,rP==rPa))
print('  full column rank over Q: %s (A\'s Lemma hypothesis)'%(rQ==n))
# how many rows are unsatisfiable no matter what (const rows)
cf=[i for i in range(len(R)) if not any(N[i]) and B[i]!=0]
print('  const-nonzero rows (always fail):',len(cf))
import sys as _s
if len(_s.argv)>2 and _s.argv[2]=='exh':
    # minimum number of rows to drop so the rest is mod-p consistent (condition (a)), exhaustive to size 6
    t0=time.time(); hit=None
    for s in range(0,7):
        cnt=0
        for D in itertools.combinations(range(len(R)),s):
            Ds=set(D)
            M2=[Naug[i] for i in range(len(R)) if i not in Ds]
            M3=[N[i] for i in range(len(R)) if i not in Ds]
            cnt+=1
            if rank_mod(M3,n,p)==rank_mod(M2,n+1,p):
                hit=(s,[R[i] for i in D]); break
        print('  cond(a) drop-size %d: %d subsets, %s (%.0fs)'%(s,cnt,'CONSISTENT mod p' if hit else 'none',time.time()-t0),flush=True)
        if hit: break
    print('  lightest mod-p-admissible violated set:',hit)
