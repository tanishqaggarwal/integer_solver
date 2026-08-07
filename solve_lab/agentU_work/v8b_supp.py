"""V8b: support propagation with the product-constraint atoms also read FORWARD
(c*(a*b) - z = 0  =>  z depends on a,b).  Same as v8 otherwise."""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); L=pickle.load(open(B+'v_leaves.pkl','rb'))
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
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
def vsof(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vsof(c,acc)
DEFSH={'(V-(V*V))','(V-(V+V))','(V-(V-V))','(V-(C*V))','(V-(V*C))','(V-(V+C))','(V-(C-V))'}
dep=collections.defaultdict(set); seeded={}
for v,lst in DEFS.items():
    for canon,sh,rv in lst:
        if sh in DEFSH: dep[find(v)] |= {find(u) for u in rv}
for v,c in CONST.items(): seeded.setdefault(find(v),set())
for sel,w,C,z,m in LEAFPIN:
    seeded[find(w)]={sel2exp[sel]}
    dep[find(z)] |= {find(sel), find(w)}
for s in sel2exp: seeded[find(s)]={sel2exp[s]}
extra=collections.Counter()
for canon,n in AT.items():
    sh=shape(n)
    tgt=None; srcs=None
    if sh in ('((V-V)-V)','((C*(V-V))-V)','((V*V)-V)','((C*(V*V))-V)'):
        tgt=n[2][1]; srcs=set(); vsof(n[1],srcs)
        if sh=='((V*V)-V)' and n[1][1][1]==n[1][2][1]: tgt=None      # boolean idempotency
    elif sh in ('((V-V)-(C*V))','((V*V)-(C*V))'):
        tgt=n[2][2][1]; srcs=set(); vsof(n[1],srcs)
    if tgt is None: continue
    ct=find(tgt)
    if ct in seeded: continue
    if ct in {find(u) for u in srcs}: continue
    dep[ct] |= {find(u) for u in srcs}; extra[sh]+=1
print('forward-read constraint atoms:',dict(extra))
for cv in seeded: dep.pop(cv,None)
cls={find(i) for i in range(NV)}
supp={c:set(seeded.get(c,())) for c in cls}
rev=collections.defaultdict(list)
for v,ds in dep.items():
    for u in ds: rev[u].append(v)
work=collections.deque(c for c in cls if supp[c]); it=0
while work:
    u=work.popleft(); it+=1
    for v in rev.get(u,()):
        if not supp[u] <= supp[v]:
            supp[v] |= supp[u]; work.append(v)
fam=collections.Counter(frozenset(supp[c]) for c in cls if supp[c])
print('propagation steps',it,' distinct nonempty supports:',len(fam))
print('by size:', sorted(collections.Counter(len(s) for s in fam).items()))
pickle.dump({'supp':supp,'par':par,'fam':dict(fam)}, open(B+'v_supp2.pkl','wb'))
