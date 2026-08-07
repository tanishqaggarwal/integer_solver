#!/usr/bin/env python3
"""AUDIT T5: are the extra knobs that F's (coarser but VERIFIED-FAITHFUL) decomposition
exposes genuine?  Test empirically: perturb each candidate and check that every equation
whose value changes lies inside R_L.  Also test exact affineness of every window atom."""
import os,sys,pickle,json,collections,random
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
from circ2 import vars_of
from fwd import compile_node
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
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<atoms>','exec')
def av(v):
    r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}}); return r
def cores(r):
    out=[0]*len(eqrows)
    for e,row in enumerate(eqrows):
        s=0
        for k,a in row:
            x=r[idx[a]]
            if x: s+=k*x
        out[e]=s
    return out
r0=av(V0); c0=cores(r0)
LEV=int(sys.argv[1]) if len(sys.argv)>1 else 0
FAIL=[12231,12270,12350,14584,18673,22044,29125]
A=set()
for e in FAIL: A|=eq_atoms[e]
for L in range(LEV):
    R=set()
    for a in A: R|=a2e[a]
    A2=set(A)
    for e in R: A2|=eq_atoms[e]
    A=A2
R=set()
for a in A: R|=a2e[a]
Vw=set()
for a in A: Vw|=avars[a]
K=sorted(u for u in Vw if v2a[u]<=A)
print('L=%d  atoms=%d  eqs=%d  vars=%d  knobs(raw,F-parse)=%d'%(LEV,len(A),len(R),len(Vw),len(K)))
bad=[]; lin_bad=[]
random.seed(11)
for u in K:
    v=list(V0); v[u]=V0[u]+1
    c1=cores(av(v))
    ch=set(e for e in range(len(eqrows)) if c1[e]!=c0[e])
    if not ch<=R: bad.append((u,sorted(ch-R)[:5]))
    # affineness: f(x+2) - f(x+1) == f(x+1) - f(x) for every equation core
    v2=list(V0); v2[u]=V0[u]+2; c2=cores(av(v2))
    for e in R:
        if c2[e]-c1[e]!=c1[e]-c0[e]: lin_bad.append((u,e)); break
print('knobs whose perturbation escapes R_L:',len(bad), bad[:5])
print('knobs on which some row of R_L is NOT affine:',len(set(x[0] for x in lin_bad)), sorted(set(x[0] for x in lin_bad))[:10])
GOOD=[u for u in K if u not in set(x[0] for x in bad)]
print('GENUINE knobs (zero collateral outside R_L):',len(GOOD))
json.dump({'L':LEV,'R':sorted(R),'K':K,'good':GOOD,'nonaffine':sorted(set(x[0] for x in lin_bad))},
          open(os.path.join(T,'window_L%d.json'%LEV),'w'))
