"""Complete fold evaluator + full-assignment constructor over the 254+1 node OR/mux tree."""
import sys, os, json, collections, pickle, re, time
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
D=pickle.load(open('ortree2.pkl','rb')); tree=dict(D['tree'])
NODE=dict(pickle.load(open('nodes.pkl','rb'))); OUT=dict(pickle.load(open('outwires.pkl','rb')))
P=pickle.load(open('pins.pkl','rb')); PIN=P['PIN']; live=P['live']; dead=P['dead']

def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen: seen.add(v); v=defrhs[v][1]
    return v
# ---- root node ----
RA,RB=D['RA'],D['RB']
prodof=collections.defaultdict(list); notof=collections.defaultdict(list); gated=collections.defaultdict(list); sumof={}
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        a,b=deref(r[1][1]),deref(r[2][1]); prodof[frozenset((a,b))].append(w)
        gated[a].append((w,b)); gated[b].append((w,a))
    elif r[0]=='-' and r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': notof[deref(r[2][1])].append(deref(w))
    elif r[0]=='+' and r[1][0]=='v' and r[2][0]=='v': sumof[frozenset((r[1][1],r[2][1]))]=w
ROOT=dict(a=RA,b=RB,
          sa=[x for nb in notof[RB] for x in prodof.get(frozenset((RA,nb)),[])],
          sb=[x for na in notof[RA] for x in prodof.get(frozenset((RB,na)),[])],
          sab=prodof.get(frozenset((RA,RB)),[]))
for k,s in (('ga','sa'),('gb','sb'),('gab','sab')): ROOT[k]=[g for x in ROOT[s] for g in gated[x]]
def outwire(N):
    res=[]
    for (pa,va) in N['ga']:
        c=[(pb,vb) for (pb,vb) in N['gb'] if frozenset((pa,pb)) in sumof]
        if len(c)!=1: return None
        pb,vb=c[0]; s2=sumof[frozenset((pa,pb))]
        c2=[(pab,vab) for (pab,vab) in N['gab'] if frozenset((s2,pab)) in sumof]
        if len(c2)!=1: return None
        pab,vab=c2[0]; res.append(dict(va=va,vb=vb,vab=vab,out=sumof[frozenset((s2,pab))]))
    return res if len(res)==2 else None
ROOTOUT=outwire(ROOT)
NODE[-1]=ROOT; OUT[-1]=ROOTOUT; tree[-1]=(RA,RB)
# ---- slot links: free wire w reading defined z ----
link={}
lr=[re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]'), re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-')]
for a in E.res:
    m=lr[0].match(a); 
    if m: u,z=int(m.group(1)),int(m.group(2))
    else:
        m=lr[1].match(a)
        if not m: continue
        u,z=int(m.group(2)),int(m.group(3))
    if u not in defrhs and z in defrhs: link[u]=z
    elif z not in defrhs and u in defrhs: link[z]=u
def chordK(A,B):
    """A,B given as (coord1,coord2); x=coord2, y=coord1"""
    ax,ay,bx,by=A[1],A[0],B[1],B[0]
    l=(by-ay)*pow(bx-ax,p-2,p)%p
    ox=(l*l-ax-bx-K)%p
    oy=(l*(ax-ox)-ay)%p
    return (oy,ox)
def build(S):
    """S: iterable of live leaves ON.  Returns (assignment dict of value wires, root value)."""
    S=set(S); v={}
    for L in S: 
        ws,Cs=PIN[L]
        v[L]=1
        for w,c in zip(ws,Cs): v[w]=c
    val={}; isl={}
    for L in tree:
        if tree[L] is None: isl[L]= (L in S); val[L]= tuple(PIN[L][1]) if L in S else None
    order=[]
    def post(n):
        if tree[n] is None: return
        for c in tree[n]: post(c)
        order.append(n)
    post(-1)
    for n in order:
        a,b=tree[n]
        la,lb=isl[a],isl[b]
        isl[n]=la or lb
        if la and lb: val[n]=chordK(val[a],val[b])
        elif la: val[n]=val[a]
        elif lb: val[n]=val[b]
        else: val[n]=None
        # assign this node's wires
        for i,d in enumerate(OUT[n]):
            for side,ch in (('va',a),('vb',b)):
                if tree[ch] is not None:      # slot wire reading child's mux output
                    v[d[side]] = (val[ch][i] if isl[ch] else 0)
            v[d['vab']] = (chordK(val[a],val[b])[i] if (la and lb) else 0)
    return v,val,isl
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    for trial,S in enumerate([[live[0]],[live[0],live[1]],rnd.sample(live,3),rnd.sample(live,17),rnd.sample(live,89),live]):
        v,val,isl=build(S)
        vv=[0]*NV
        for k,x in v.items(): vv[k]=x
        r=E.run(vv)
        nz=[i for i,x in enumerate(r) if x%p]
        print('|S|=%-4d  nonzero residual atoms mod p: %d'%(len(S),len(nz)))
        if nz and len(nz)<=8:
            for i in nz: print('     ',E.res[i][:150])
        print('     root value =',val[-1])
