import sys, os, json, re, collections, pickle
F='/home/user/integer_solver/solve_lab/agentV_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
# uses: var -> list of defined vars whose rhs mentions it
uses=collections.defaultdict(list)
for w,r in defrhs.items():
    for u in vars_of(r): uses[u].append(w)

def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen:
        seen.add(v); v=defrhs[v][1]
    return v

def as_or(v):
    v=deref(v); r=defrhs.get(v)
    if r is None or r[0]!='-': return None
    def unv(n):
        if n[0]=='v':
            w=deref(n[1]); return defrhs.get(w)
        return n
    Ln=unv(r[1]); Rn=unv(r[2])
    if Ln is None or Rn is None: return None
    if Ln[0]=='+' and Rn[0]=='*':
        s1={node_str(Ln[1]),node_str(Ln[2])}; s2={node_str(Rn[1]),node_str(Rn[2])}
        if s1==s2: return (Ln[1],Ln[2])
    return None

def as_not(v):
    """x = 1 - y"""
    v=deref(v); r=defrhs.get(v)
    if r is None or r[0]!='-': return None
    if r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': return deref(r[2][1])
    return None

def as_mul(v):
    v=deref(v); r=defrhs.get(v)
    if r is None or r[0]!='*': return None
    a,b=r[1],r[2]
    if a[0]=='v' and b[0]=='v': return (deref(a[1]),deref(b[1]))
    return None

# --- build tree ---
tree={}   # node -> (child_a, child_b) in OR-var terms, or None for leaf
def build(v):
    v=deref(v)
    if v in tree: return v
    o=as_or(v)
    if o is None:
        tree[v]=None; return v
    ca=build(o[0][1]) if o[0][0]=='v' else None
    cb=build(o[1][1]) if o[1][0]=='v' else None
    tree[v]=(ca,cb); return v

if __name__=='__main__':
    R_A=build(8599); R_B=build(21839)
    print('OR nodes total', sum(1 for v in tree if tree[v]), 'leaves', sum(1 for v in tree if tree[v] is None))
    # find selectors: for each OR node with children a,b, look for products
    # need the OR-var of each child; leaf's "live" indicator is the leaf var itself
    sel={}
    notof={}   # y -> x with x = 1-y
    for w in defrhs:
        n=as_not(w)
        if n is not None: notof.setdefault(n,[]).append(deref(w))
    prodof=collections.defaultdict(list)
    for w in defrhs:
        m=as_mul(w)
        if m is not None: prodof[frozenset(m)].append(deref(w))
    cnt=collections.Counter()
    selmap={}
    for n,ch in tree.items():
        if ch is None: continue
        a,b=ch
        na=notof.get(a,[]); nb=notof.get(b,[])
        s_ab=prodof.get(frozenset((a,b)),[])
        s_a = [x for nbv in nb for x in prodof.get(frozenset((a,nbv)),[])]
        s_b = [x for nav in na for x in prodof.get(frozenset((b,nav)),[])]
        selmap[n]=dict(a=a,b=b,s_a=s_a,s_b=s_b,s_ab=s_ab)
        cnt[(len(s_a)>0,len(s_b)>0,len(s_ab)>0)]+=1
    print('selector presence (a-only,b-only,both):',cnt.most_common())
    pickle.dump({'tree':tree,'selmap':selmap,'RA':R_A,'RB':R_B},open('ortree2.pkl','wb'))
    # depth profile
    def depth(v):
        if tree[v] is None: return 0
        return 1+max(depth(tree[v][0]),depth(tree[v][1]))
    print('depth A',depth(R_A),'depth B',depth(R_B))
