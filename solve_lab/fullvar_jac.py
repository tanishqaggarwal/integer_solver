#!/usr/bin/env python3
"""DECISIVE: mod-p feasibility of repairing the 11 fails using ALL 38748 variables
(off-manifold), instead of only the 8583 forward-eval free inputs.
Prior analyses restricted to free inputs and found a dim-1 obstruction (defect 1).
Off-manifold, gate outputs are independent variables -> far more freedom."""
import sys, os, re, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
p=H.p
NV=38748
class D:
    __slots__=('v','g')
    def __init__(s,v,g=None): s.v=v%p; s.g=g if g is not None else {}
    def __add__(s,o):
        if isinstance(o,D):
            g=dict(s.g)
            for k,c in o.g.items():
                n=(g.get(k,0)+c)%p
                if n: g[k]=n
                elif k in g: del g[k]
            return D(s.v+o.v,g)
        return D(s.v+o,dict(s.g))
    __radd__=__add__
    def __sub__(s,o):
        if isinstance(o,D):
            g=dict(s.g)
            for k,c in o.g.items():
                n=(g.get(k,0)-c)%p
                if n: g[k]=n
                elif k in g: del g[k]
            return D(s.v-o.v,g)
        return D(s.v-o,dict(s.g))
    def __rsub__(s,o): return D(o-s.v,{k:(-c)%p for k,c in s.g.items()})
    def __mul__(s,o):
        if isinstance(o,D):
            g={}
            for k,c in s.g.items():
                n=(c*o.v)%p
                if n: g[k]=n
            for k,c in o.g.items():
                n=(g.get(k,0)+c*s.v)%p
                if n: g[k]=n
                elif k in g: del g[k]
            return D(s.v*o.v,g)
        g={}
        for k,c in s.g.items():
            n=(c*o)%p
            if n: g[k]=n
        return D(s.v*o,g)
    __rmul__=__mul__
    def __neg__(s): return D(-s.v,{k:(-c)%p for k,c in s.g.items()})
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
base=[H.val[i] for i in range(NV)]
pat=re.compile(r'x_(\d+)')
print('compiling equations...',flush=True); t0=time.time()
comp=[]
for L in open('/home/user/integer_solver/EQUATIONS.txt'):
    L=L.strip()
    if not L: continue
    if L.endswith('= 0'): L=L[:-3].rstrip()
    comp.append(compile(pat.sub(lambda m:'v['+m.group(1)+']',L),'<e>','eval'))
print(f'  {len(comp)} eqs, {time.time()-t0:.0f}s',flush=True)
vd=[D(base[i],{i:1}) for i in range(NV)]
ns={'v':vd,'__builtins__':{}}
print('evaluating jacobian rows (all 38748 columns)...',flush=True); t0=time.time()
rows=[]; nfail=0
for i,c in enumerate(comp):
    r=eval(c,ns)
    if isinstance(r,D): val,g=r.v,r.g
    else: val,g=r%p,{}
    if val!=0: nfail+=1
    if g or val!=0: rows.append((g,(-val)%p))
    if i%8000==0: print(f'   eq {i}, {time.time()-t0:.0f}s',flush=True)
print(f'  done {time.time()-t0:.0f}s; constraining rows={len(rows)}; nonzero residual rows={nfail}',flush=True)
def inv(a): return pow(a%p,p-2,p)
print('sparse GE (Rouche-Capelli over ALL variables)...',flush=True); t0=time.time()
piv={}; rank=0; incons=0
for rd0,rb0 in rows:
    rd=dict(rd0); rb=rb0
    stack=list(rd)
    while stack:
        c=stack.pop()
        if c not in rd or c not in piv: continue
        f=rd[c]
        if f==0: del rd[c]; continue
        prd,prb=piv[c]
        for k,val2 in prd.items():
            n=(rd.get(k,0)-f*val2)%p
            if n:
                if k not in rd: stack.append(k)
                rd[k]=n
            elif k in rd: del rd[k]
        rb=(rb-f*prb)%p
    rd={k:v for k,v in rd.items() if v}
    if not rd:
        if rb%p!=0: incons+=1
        continue
    c=next(iter(rd)); ic=inv(rd[c])
    piv[c]=({k:(v*ic)%p for k,v in rd.items()},(rb*ic)%p); rank+=1
    if rank%4000==0: print(f'   rank={rank}, {time.time()-t0:.0f}s',flush=True)
print(f'\nRESULT: rank={rank}, inconsistent_rows={incons}, {time.time()-t0:.0f}s',flush=True)
if incons==0:
    print('*** CONSISTENT over all 38748 variables -> an off-manifold mod-p repair EXISTS ***',flush=True)
else:
    print(f'*** INCONSISTENT ({incons} rows) even off-manifold ***',flush=True)
