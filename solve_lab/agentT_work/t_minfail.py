#!/usr/bin/env python3
"""AUDIT T6: exact minfail over the ENLARGED knob set that F's faithful (coarser)
decomposition exposes.  Agent A's Theorem A/B is proved over 9 knobs at L=0 / 109 at L=6;
this window has 24 / 235.  If minfail < 7 here, 39,026 is beatable."""
import os,sys,pickle,json,collections,itertools,time
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
I=os.path.abspath(os.path.join(T,'..','agentI_work')); sys.path.insert(0,I)
from circ2 import vars_of
from fwd import compile_node
from intsolve import solve_int
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
def cores(v):
    r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}})
    out={}
    for e,row in enumerate(eqrows):
        s=0
        for k,a in row:
            x=r[idx[a]]
            if x: s+=k*x
        out[e]=s
    return out
LEV=int(sys.argv[1]); BUD=int(sys.argv[2]) if len(sys.argv)>2 else 6
FAIL=[12231,12270,12350,14584,18673,22044,29125]
A=set()
for e in FAIL: A|=eq_atoms[e]
for L in range(LEV):
    R=set()
    for a in A: R|=a2e[a]
    A2=set(A)
    for e in R: A2|=eq_atoms[e]
    A=A2
R=sorted({e for a in A for e in a2e[a]})
Vw=set()
for a in A: Vw|=avars[a]
K=sorted(u for u in Vw if v2a[u]<=A)
print('L=%d atoms=%d eqs=%d knobs=%d'%(LEV,len(A),len(R),len(K)),flush=True)
t0=time.time(); c0=cores(V0)
Mat={e:[] for e in R}
for u in K:
    v=list(V0); v[u]=V0[u]+1; c1=cores(v)
    for e in R: Mat[e].append(c1[e]-c0[e])
print('jacobian built %.0fs'%(time.time()-t0),flush=True)
rows=[];bvec=[];triv_fail=[];triv_ok=[]
for e in R:
    if any(Mat[e]):
        rows.append((e,Mat[e],-c0[e]))
    elif c0[e]!=0: triv_fail.append(e)
    else: triv_ok.append(e)
print('nontrivial rows=%d  forced-fail(const!=0)=%d  const-ok=%d'%(len(rows),len(triv_fail),len(triv_ok)),flush=True)
print('rows currently violated at S:',sum(1 for e,_,b in rows if b!=0)+len(triv_fail),flush=True)
n=len(rows); best=None
for s in range(0,BUD+1):
    if len(triv_fail)+s>BUD: break
    found=False
    cnt=0
    for D in itertools.combinations(range(n),s):
        Ds=set(D)
        Am=[rows[i][1] for i in range(n) if i not in Ds]
        bm=[rows[i][2] for i in range(n) if i not in Ds]
        cnt+=1
        if solve_int(Am,bm) is not None:
            best=len(triv_fail)+s
            print('SOLVABLE sacrificing %d rows: %s  => total failing %d'%(s,[rows[i][0] for i in D],best),flush=True)
            found=True; break
    print('  size %d: %d subsets tested %s (%.0fs)'%(s,cnt,'HIT' if found else 'none',time.time()-t0),flush=True)
    if found: break
print('minfail over this window =',best if best is not None else '> %d'%(len(triv_fail)+BUD),flush=True)
