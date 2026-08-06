#!/usr/bin/env python3
"""The 11 failing equations involve only 75 variables; exactly 1099 equations touch them.
Modifying ONLY those 75 leaves all other 37934 equations untouched -> solving these 1099
gives a COMPLETE solution. Test mod-p consistency of that local system (off-manifold:
all 75 are independent unknowns, no forward-eval)."""
import sys, os, re, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
os.chdir('/home/user/integer_solver/solve_lab')
import heal_harness as H
p=H.p; NV=38748
lines=open('/home/user/integer_solver/EQUATIONS.txt').read().splitlines()
pat=re.compile(r'x_(\d+)')
eqvars=[set(int(m) for m in pat.findall(L)) for L in lines]
var_eqs=collections.defaultdict(set)
for i,vs in enumerate(eqvars):
    for v in vs: var_eqs[v].add(i)
fails=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
import sys as _s
RING=int(_s.argv[1]) if len(_s.argv)>1 else 1
V=set()
for e in fails: V|=eqvars[e]
E=set().union(*[var_eqs[v] for v in V])
for _ in range(RING-1):
    for e in list(E): V|=eqvars[e]
    E=set().union(*[var_eqs[v] for v in V])
V=sorted(V); E=sorted(E)
print(f'local system: {len(V)} variables, {len(E)} equations')
col={v:k for k,v in enumerate(V)}
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
base=[H.val[i] for i in range(NV)]
comp={i:compile(pat.sub(lambda m:'v['+m.group(1)+']',(lines[i][:-3].rstrip() if lines[i].strip().endswith('= 0') else lines[i])),'<e>','eval') for i in E}
class D:
    __slots__=('v','g')
    def __init__(s,v,g=None): s.v=v%p; s.g=g or {}
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
def inv(a): return pow(a%p,p-2,p)
def build():
    vd=[None]*NV
    for i in range(NV): vd[i]=base[i]
    for v in V: vd[v]=D(base[v],{col[v]:1})
    ns={'v':vd,'__builtins__':{}}
    rows=[]
    for i in E:
        r=eval(comp[i],ns)
        if isinstance(r,D): val,g=r.v,r.g
        else: val,g=r%p,{}
        if g or val: rows.append((g,(-val)%p,i))
    return rows
rows=build()
nz=[i for g,b,i in rows if b!=0]
print(f'constraining rows={len(rows)}, rows with nonzero residual={len(nz)} -> {nz}')
piv={}; rank=0; incons=[]
for rd0,rb0,ei in rows:
    rd=dict(rd0); rb=rb0
    while True:
        P=[c for c in rd if c in piv]
        if not P: break
        c=min(P); f=rd[c]
        prd,prb=piv[c]
        for k,val2 in prd.items():
            n=(rd.get(k,0)-f*val2)%p
            if n: rd[k]=n
            elif k in rd: del rd[k]
        rb=(rb-f*prb)%p
    rd={k:v for k,v in rd.items() if v}
    if not rd:
        if rb%p: incons.append(ei)
        continue
    c=min(rd); ic=inv(rd[c])
    piv[c]=({k:(v*ic)%p for k,v in rd.items()},(rb*ic)%p); rank+=1
print(f'\nRANK={rank} of {len(V)} columns; INCONSISTENT rows={len(incons)}')
if not incons:
    print('*** LOCAL SYSTEM CONSISTENT mod p -> off-manifold repair exists in 75 vars ***')
else:
    print('*** inconsistent at:',incons[:10])
