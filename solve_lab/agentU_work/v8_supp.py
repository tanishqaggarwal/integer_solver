"""U8: forward selector-support on the canonical definition DAG (my parse only)."""
import pickle, collections, sys
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb'))
L=pickle.load(open(B+'v_leaves.pkl','rb'))
AT=pickle.load(open(B+'v_atoms.pkl','rb'))['AT']
DEFS=D['DEFS']; COPY=D['COPY']; CONST=D['CONST']; LEAFPIN=D['LEAFPIN']
sel2exp=L['sel2exp']; NV=38748
par=list(range(NV))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def uni(a,b):
    a,b=find(a),find(b)
    if a!=b: par[a]=b
for a,b in COPY: uni(a,b)
DEFSH={'(V-(V*V))','(V-(V+V))','(V-(V-V))','(V-(C*V))','(V-(V*C))','(V-(V+C))','(V-(C-V))'}
dep=collections.defaultdict(set)     # canonical v -> canonical rhs deps
seeded={}
for v,lst in DEFS.items():
    cv=find(v)
    for canon,sh,rv in lst:
        if sh in ('(V-C)',): seeded.setdefault(cv,set())
        elif sh in DEFSH: dep[cv] |= {find(u) for u in rv}
for v,c in CONST.items(): seeded.setdefault(find(v),set())
# leaf pins: sel*(w-C) - m*z   -> w carries the leaf, z is gated
for sel,w,C,z,m in LEAFPIN:
    seeded[find(w)]={sel2exp[sel]}
    dep[find(z)] |= {find(sel), find(w)}
for s in sel2exp: seeded[find(s)]={sel2exp[s]}
dep.pop(0,None)
for cv in seeded: dep.pop(cv,None)      # seeds override
# also treat  ((V-V)-V) and ((C*(V-V))-V) and ((V-V)-(C*V)) as defs of the trailing var
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
def vsof(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vsof(c,acc)
extra=0
for canon,n in AT.items():
    sh=shape(n)
    if sh in ('((V-V)-V)','((C*(V-V))-V)','((V-V)-(C*V))'):
        tgt = n[2][1] if n[2][0]=='var' else n[2][2][1]
        lhs=set(); vsof(n[1],lhs)
        ct=find(tgt)
        if ct in seeded: continue
        dep[ct] |= {find(u) for u in lhs}; extra+=1
print('extra defs from residual-shaped atoms:',extra)
cls={find(i) for i in range(NV)}
print('classes',len(cls),'seeded',len(seeded),'with deps',len(dep),'undefined',len(cls-set(seeded)-set(dep)))
# iterative fixpoint (handles cycles: least fixpoint = union closure, monotone)
supp={c:set(seeded.get(c,())) for c in cls}
work=collections.deque(c for c in cls if supp[c])
rev=collections.defaultdict(list)
for v,ds in dep.items():
    for u in ds: rev[u].append(v)
it=0
while work:
    u=work.popleft(); it+=1
    for v in rev.get(u,()):
        if not supp[u] <= supp[v]:
            supp[v] |= supp[u]; work.append(v)
print('propagation steps',it)
h=collections.Counter(len(supp[c]) for c in cls)
print('support-size histogram:', sorted(h.items())[:6],'...',sorted(h.items())[-8:])
fam=collections.Counter(frozenset(supp[c]) for c in cls if 0<len(supp[c]))
print('distinct nonempty supports:', len(fam))
sz=collections.Counter(len(s) for s in fam)
print('by size:', sorted(sz.items()))
pickle.dump({'supp':supp,'par':par,'fam':{s:n for s,n in fam.items()}}, open(B+'v_supp.pkl','wb'))
