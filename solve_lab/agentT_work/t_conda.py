#!/usr/bin/env python3
"""AUDIT T7: agent A's condition (a) -- minimum-weight mod-p-admissible violated set --
recomputed over the ENLARGED knob set (24/88/166/235 knobs instead of A's 9/32/80/109).
A's claim: 'the lightest weight the filter admits is 7 at every depth'."""
import os,sys,pickle,json,collections,itertools,time
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
LEV=int(sys.argv[1]); MAXW=int(sys.argv[2]) if len(sys.argv)>2 else 6
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
m=len(R); n=len(K)
N=[[0]*n for _ in range(m)]
for j,u in enumerate(K):
    v=list(V0); v[u]=V0[u]+1; c1=cores(v)
    for i in range(m): N[i][j]=c1[i]-c0[i]
B=[(-c0[i])%p for i in range(m)]
# left kernel of N mod p: rows lam with lam*N = 0
Aug=[[N[i][j]%p for j in range(n)]+[1 if k==i else 0 for k in range(m)] for i in range(m)]
piv=0; pivcols=[]
for c in range(n):
    k=None
    for i in range(piv,m):
        if Aug[i][c]: k=i;break
    if k is None: continue
    Aug[piv],Aug[k]=Aug[k],Aug[piv]
    inv=pow(Aug[piv][c],p-2,p); Aug[piv]=[(x*inv)%p for x in Aug[piv]]
    for i in range(m):
        if i!=piv and Aug[i][c]:
            f=Aug[i][c]; Aug[i]=[(Aug[i][t]-f*Aug[piv][t])%p for t in range(n+m)]
    pivcols.append(c); piv+=1
Wb=[Aug[i][n:] for i in range(piv,m)]      # w x m
w=len(Wb)
g=[sum(Wb[r][i]*B[i] for i in range(m))%p for r in range(w)]
print('L=%d rows=%d knobs=%d rank_p(N)=%d  left-kernel dim w=%d'%(LEV,m,n,piv,w),flush=True)
print('syndrome g == 0 ? %s  (if 0, ALL rows are already mod-p consistent)'%all(x==0 for x in g),flush=True)
cols=[[Wb[r][i] for r in range(w)] for i in range(m)]      # column i of Wb
# min |D| with g in span{cols[i] : i in D}: DFS over increasing size
t0=time.time(); best=None
def reduce_vec(v,basis):
    v=v[:]
    for (lead,bv) in basis:
        if v[lead]:
            f=v[lead]
            v=[(v[t]-f*bv[t])%p for t in range(w)]
    return v
def add_basis(v,basis):
    v=reduce_vec(v,basis)
    lead=next((t for t in range(w) if v[t]),None)
    if lead is None: return None
    inv=pow(v[lead],p-2,p); v=[(x*inv)%p for x in v]
    return basis+[(lead,v)]
nodes=[0]
def dfs(start,basis,depth,chosen,budget):
    global best
    nodes[0]+=1
    if all(x==0 for x in reduce_vec(g,basis)):
        best=list(chosen); return True
    if depth==budget: return False
    for i in range(start,m):
        if m-i < budget-depth: break
        nb=add_basis(cols[i],basis)
        if nb is None: nb=basis
        chosen.append(i)
        if dfs(i+1,nb,depth+1,chosen,budget): return True
        chosen.pop()
    return False
for budget in range(0,MAXW+1):
    nodes[0]=0
    if dfs(0,[],0,[],budget):
        print('  weight %d: ADMISSIBLE  D=%s  (%d nodes, %.0fs)'%(budget,[R[i] for i in best],nodes[0],time.time()-t0),flush=True)
        break
    print('  weight %d: none admissible (%d nodes, %.0fs)'%(budget,nodes[0],time.time()-t0),flush=True)
else:
    print('  no admissible set of weight <= %d'%MAXW,flush=True)
