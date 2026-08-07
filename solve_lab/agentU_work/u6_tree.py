"""U6: recover the liveness (OR) tree over the 256 leaf selectors, from my own parse only."""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'u_defs.pkl','rb'))
L=pickle.load(open(B+'u_leaves.pkl','rb'))
DEFS=D['DEFS']; COPY=D['COPY']; CONST=D['CONST']
sel2exp=L['sel2exp']
NV=38748
par=list(range(NV))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def uni(a,b):
    a,b=find(a),find(b)
    if a!=b: par[a]=b
for a,b in COPY: uni(a,b)
print('copy classes:', len({find(i) for i in range(NV)}))
# canonical defs
CD=collections.defaultdict(list)
for v,lst in DEFS.items():
    for canon,sh,rv in lst:
        CD[find(v)].append((sh,tuple(sorted(find(u) for u in rv)),canon,v))
# index: rhs pattern -> var
ADD={}; MUL={}; SUB=[]
for v,lst in CD.items():
    for sh,rv,canon,orig in lst:
        if sh=='(V-(V+V))' and len(rv)==2: ADD.setdefault(rv,set()).add(v)
        if sh=='(V-(V*V))' and len(rv)==2: MUL.setdefault(rv,set()).add(v)
        if sh=='(V-(V-V))': SUB.append((v,canon,orig))
# also parse the exact operand order for the (V-(V-V)) shape
import re
ATOMS=pickle.load(open(B+'u_atoms.pkl','rb'))['ATOMS']
subops={}
for v,canon,orig in SUB:
    n=ATOMS[canon]           # ('sub', var, ('sub', X, Y))
    subops.setdefault(v,[]).append((find(n[2][1][1]), find(n[2][2][1])))
LIVE={find(s):('leaf',s) for s in sel2exp}
print('leaf live classes',len(LIVE))
changed=True; internal={}
while changed:
    changed=False
    for v,ops in subops.items():
        if v in LIVE or v in internal: continue
        for (t1,t2) in ops:
            # need t1 = a+b, t2 = a*b, a,b live
            a1=[k for k in ADD if t1 in ADD[k]]
            a2=[k for k in MUL if t2 in MUL[k]]
            hit=[k for k in a1 if k in a2 and len(k)==2 and k[0] in LIVE and k[1] in LIVE and k[0]!=k[1]]
            if hit:
                internal[v]=hit[0]; LIVE[v]=('or',hit[0]); changed=True; break
print('internal OR nodes found:', len(internal))
# support closure
supp={}
def get(v):
    if v in supp: return supp[v]
    kind=LIVE[v]
    if kind[0]=='leaf': r=frozenset([sel2exp[kind[1]]])
    else:
        a,b=kind[1]; r=get(a)|get(b)
    supp[v]=r; return r
roots=[v for v in LIVE if all(v not in internal[w] for w in internal)]
for v in LIVE: get(v)
sizes=collections.Counter(len(supp[v]) for v in LIVE)
print('support size histogram (top):', sizes.most_common(12))
print('max support', max(len(s) for s in supp.values()))
# who is the root
tops=[v for v in LIVE if len(supp[v])==max(len(s) for s in supp.values())]
print('top nodes', tops, [len(supp[t]) for t in tops])
pickle.dump({'LIVE':LIVE,'internal':internal,'supp':{k:set(v) for k,v in supp.items()},'find':par}, open(B+'u_tree.pkl','wb'))
