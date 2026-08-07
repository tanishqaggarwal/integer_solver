import pickle, re, ast, collections, sys
m=pickle.load(open('model3.pkl','rb'))
atoms=m['atoms']; eqt=m['eq_terms']
VAR=re.compile(r'x_(\d+)')
def vars_of(s): return set(int(x) for x in VAR.findall(s))
# improved def detection: solve atom for a variable appearing linearly with coeff +-1
# forms:  X - RHS ;  RHS - X ; X - A - B
info=[]
for i,a in enumerate(atoms):
    t=ast.parse(a,mode='eval').body
    d=None
    # peel  A - B  chains: left-assoc
    if isinstance(t,ast.BinOp) and isinstance(t.op,ast.Sub):
        # collect  t = L - R
        L,R=t.left,t.right
        if isinstance(L,ast.Name):
            ov=int(L.id[2:]); rv=vars_of(ast.unparse(R))
            if ov not in rv: d=(ov,'L')
        elif isinstance(R,ast.Name):
            ov=int(R.id[2:]); lv=vars_of(ast.unparse(L))
            if ov not in lv: d=(ov,'R')
        if d is None and isinstance(L,ast.BinOp) and isinstance(L.op,ast.Sub) and isinstance(L.left,ast.Name):
            ov=int(L.left.id[2:]); rv=vars_of(ast.unparse(L.right))|vars_of(ast.unparse(R))
            if ov not in rv: d=(ov,'LL')
    elif isinstance(t,ast.Name):
        d=(int(t.id[2:]),'bare')
    info.append((d, vars_of(a)))
ndef=sum(1 for d,_ in info if d)
print("defs",ndef,"of",len(atoms))
defs=collections.defaultdict(list)
for i,(d,vs) in enumerate(info):
    if d: defs[d[0]].append(i)
print("vars defined",len(defs))
print("hist", sorted(collections.Counter(len(v) for v in defs.values()).items()))
allv=set()
for _,vs in info: allv|=vs
print("vars appearing in atoms:",len(allv))
free=allv-set(defs)
print("free (never defined):",len(free))
# check acyclicity using ONE def per var (first)
import sys
sys.setrecursionlimit(100000)
defatom={v:ids[0] for v,ids in defs.items()}
color={}
cyc=[]
def dfs(v):
    st=[(v,0)]
    while st:
        u,ph=st.pop()
        if ph==0:
            c=color.get(u,0)
            if c==1: cyc.append(u); continue
            if c==2: continue
            color[u]=1
            st.append((u,1))
            ai=defatom.get(u)
            if ai is not None:
                for w in info[ai][1]:
                    if w!=u: st.append((w,0))
        else:
            color[u]=2
for v in list(defs): 
    if color.get(v,0)==0: dfs(v)
print("cycle hits:",len(cyc), cyc[:20])
pickle.dump({'info':info,'defs':dict(defs),'free':free}, open('dag.pkl','wb'))
