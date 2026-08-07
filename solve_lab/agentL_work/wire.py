import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
uses=collections.defaultdict(list)
for w,r in defrhs.items():
    for u in vars_of(r): uses[u].append(w)
D=pickle.load(open('ortree2.pkl','rb'))
tree=D['tree']; selmap=D['selmap']
NODE=pickle.load(open('nodes.pkl','rb'))

# find the sum wire that combines two given product wires
sumof={}
for w,r in defrhs.items():
    if r[0]=='+' and r[1][0]=='v' and r[2][0]=='v':
        sumof[frozenset((r[1][1],r[2][1]))]=w

def outwire(n):
    """Return dict coord-> (prod_a, prod_b, prod_ab, sum2, out) for the node's mux, matching coords."""
    N=NODE[n]
    res=[]
    for (pa,va) in N['ga']:
        # find pb sharing the sum with pa
        cand=[(pb,vb) for (pb,vb) in N['gb'] if frozenset((pa,pb)) in sumof]
        if len(cand)!=1: return None
        pb,vb=cand[0]; s2=sumof[frozenset((pa,pb))]
        cand2=[(pab,vab) for (pab,vab) in N['gab'] if frozenset((s2,pab)) in sumof]
        if len(cand2)!=1: return None
        pab,vab=cand2[0]; out=sumof[frozenset((s2,pab))]
        res.append(dict(va=va,vb=vb,vab=vab,out=out))
    return res if len(res)==2 else None

ok=0; bad=[]
OUT={}
for n in NODE:
    r=outwire(n)
    if r is None: bad.append(n)
    else: ok+=1; OUT[n]=r
print('nodes with clean 2-coord 3-way mux:',ok,'  failed:',len(bad),bad[:5])

# consistency: for node n with OR-child c, does n read c's output wires?
mism=0; match=0; leafread=0
for n,N in NODE.items():
    if n not in OUT: continue
    for side,ch in (('ga',N['a']),('gb',N['b'])):
        vals={d['va'] if side=='ga' else d['vb'] for d in OUT[n]}
        if tree[ch] is not None:
            couts={d['out'] for d in OUT.get(ch,[])}
            if couts==vals: match+=1
            else: mism+=1; 
        else: leafread+=1
print('OR-child value wires == child mux output: match',match,'mismatch',mism,' leaf-child reads',leafread)
pickle.dump(OUT,open('outwires.pkl','wb'))
