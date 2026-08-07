import sys, os, json, re, collections
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)

def deref(v):
    """follow chains x := y"""
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen:
        seen.add(v); v=defrhs[v][1]
    return v

def as_or(v):
    """if x_v = a+b-a*b return (a,b) else None"""
    v=deref(v)
    r=defrhs.get(v)
    if r is None or r[0]!='-': return None
    L,R=r[1],r[2]
    def unv(n):
        if n[0]=='v': 
            w=deref(n[1]); return defrhs.get(w), w
        return n,None
    Ln,Lw=unv(L); Rn,Rw=unv(R)
    if Ln is None or Rn is None: return None
    if Ln[0]=='+' and Rn[0]=='*':
        a1,b1=Ln[1],Ln[2]; a2,b2=Rn[1],Rn[2]
        s1={node_str(a1),node_str(b1)}; s2={node_str(a2),node_str(b2)}
        if s1==s2:
            return (a1,b1)
    return None

def leaves(v, out, ors, depth=0):
    o=as_or(v)
    if o is None:
        out.append(deref(v)); return
    ors.append((deref(v),depth))
    for c in o:
        if c[0]=='v': leaves(c[1],out,ors,depth+1)
        else: out.append(('lit',node_str(c)))

if __name__=='__main__':
    root=int(sys.argv[1])
    out=[];ors=[]
    leaves(root,out,ors)
    print('root x%d: %d OR-leaves, %d OR-nodes'%(root,len(out),len(ors)))
    print('depth hist', sorted(collections.Counter(d for _,d in ors).items()))
    free=[v for v in out if v not in defrhs]
    print('leaves that are FREE vars: %d / %d'%(len(free),len(out)))
    print('first leaves', out[:12])
    for v in out[:5]:
        print(' x%d res-atoms: %s'%(v,[a[:100] for a in resby.get(v,[])]))
    json.dump({'root':root,'leaves':[int(x) for x in out],'ornodes':[[int(a),d] for a,d in ors]},open('ortree_%d.json'%root,'w'))
