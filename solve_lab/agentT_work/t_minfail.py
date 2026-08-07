#!/usr/bin/env python3
"""AUDIT T6: exact minfail over the ENLARGED knob set exposed by F's faithful (coarser)
decomposition.  Agent A's Theorem A/B is proved over 9 knobs at L=0, 32 at L=2, 109 at L=6.
The same windows in F's parse carry 24 / 88 / 235 knobs, all verified zero-collateral and
exactly affine (t_knobs.py).  Re-run the optimisation with the exclusion removed."""
import os,sys,pickle,json,collections,time
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
I=os.path.abspath(os.path.join(T,'..','agentI_work')); sys.path.insert(0,I)
from circ2 import vars_of
from fwd import compile_node
from eq8680 import minfail_bnb
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
LEV=int(sys.argv[1]); BUD=int(sys.argv[2]) if len(sys.argv)>2 else 6
TL=int(sys.argv[3]) if len(sys.argv)>3 else 3000
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
print('L=%d atoms=%d eqs=%d vars=%d knobs=%d'%(LEV,len(A),len(R),len(Vw),len(K)),flush=True)
# only atoms of the window matter for rows of R? no: rows of R contain foreign atoms too.
need=sorted({a for e in R for a in eq_atoms[e]})
sub=[names[i] for i in need]; pos={a:i for i,a in enumerate(need)}
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in sub)+']','<at>','exec')
rr=[0]*len(need)
def cores(v):
    exec(prog,{'v':v,'r':rr,'__builtins__':{}})
    out=[]
    for e in R:
        s=0
        for k,a in eqrows[e]:
            x=rr[pos[idx[a]]]
            if x: s+=k*x
        out.append(s)
    return out
t0=time.time(); c0=cores(V0)
print('rows violated at S: %d'%sum(1 for x in c0 if x),flush=True)
Mat=[[] for _ in R]
for u in K:
    v=list(V0); v[u]=V0[u]+1; c1=cores(v)
    for i in range(len(R)): Mat[i].append(c1[i]-c0[i])
print('jacobian %d x %d built in %.0fs'%(len(R),len(K),time.time()-t0),flush=True)
mf,forced,nact,nodes=minfail_bnb(list(range(len(R))),c0,Mat,budget=BUD,tlimit=TL)
print('forced(const-nonzero rows)=%d  active=%d  nodes=%d'%(forced,nact,nodes),flush=True)
if mf is None: print('RESULT: minfail > %d  -> cannot beat 39,026 in this window'%BUD,flush=True)
elif mf=='timeout': print('RESULT: TIMEOUT after %d nodes'%nodes,flush=True)
else: print('RESULT: minfail = %d %s'%(mf,'*** BEATS 39,026 ***' if mf<7 else ''),flush=True)
