"""Complete tree model: 256 leaves (178 live w/ constants, 78 dead), 254 OR nodes + root.
For every node determine: children, selectors, per-coord (va,vb,vab,out), and slot links."""
import sys, os, json, collections, pickle, re
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
D=pickle.load(open('ortree2.pkl','rb')); tree=D['tree']
NODE=pickle.load(open('nodes.pkl','rb')); OUT=pickle.load(open('outwires.pkl','rb'))

def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen: seen.add(v); v=defrhs[v][1]
    return v

# ---- leaves ----
leaves=[v for v in tree if tree[v] is None]
live={}; dead=[]
pinre=re.compile(r'^\(\(x(\d+)\*\(x(\d+)-(\d+)\)\)-')
for L in leaves:
    if L in defrhs: dead.append(L); continue
    pins={}
    for a in resby[L]:
        m=pinre.match(a)
        if m and int(m.group(1))==L: pins[int(m.group(2))]=int(m.group(3))%p
    assert len(pins)==2,(L,resby[L])
    live[L]=pins
print('leaves %d  live %d  dead %d'%(len(leaves),len(live),len(dead)))
assert all(defrhs[d]==('c',0) for d in dead), 'dead leaves not literal 0'
print('all %d dead leaves are literally := 0  OK'%len(dead))

# ---- root node: find the OR node above x8599 / x21839 ----
RA,RB=8599,21839
# find products sel = a*b etc at root
prodof=collections.defaultdict(list)
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        prodof[frozenset((deref(r[1][1]),deref(r[2][1])))].append(w)
notof=collections.defaultdict(list)
for w,r in defrhs.items():
    if r[0]=='-' and r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': notof[deref(r[2][1])].append(deref(w))
s_ab=prodof.get(frozenset((RA,RB)),[])
s_a=[x for nb in notof[RB] for x in prodof.get(frozenset((RA,nb)),[])]
s_b=[x for na in notof[RA] for x in prodof.get(frozenset((RB,na)),[])]
print('ROOT selectors  a-only',s_a,' b-only',s_b,' both',s_ab)
gated=collections.defaultdict(list)
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        gated[deref(r[1][1])].append((w,deref(r[2][1]))); gated[deref(r[2][1])].append((w,deref(r[1][1])))
sumof={}
for w,r in defrhs.items():
    if r[0]=='+' and r[1][0]=='v' and r[2][0]=='v': sumof[frozenset((r[1][1],r[2][1]))]=w
ROOT=dict(a=RA,b=RB,sa=s_a,sb=s_b,sab=s_ab,
          ga=[g for s in s_a for g in gated[s]],gb=[g for s in s_b for g in gated[s]],
          gab=[g for s in s_ab for g in gated[s]])
def outwire(N):
    res=[]
    for (pa,va) in N['ga']:
        c=[(pb,vb) for (pb,vb) in N['gb'] if frozenset((pa,pb)) in sumof]
        if len(c)!=1: return None
        pb,vb=c[0]; s2=sumof[frozenset((pa,pb))]
        c2=[(pab,vab) for (pab,vab) in N['gab'] if frozenset((s2,pab)) in sumof]
        if len(c2)!=1: return None
        pab,vab=c2[0]
        res.append(dict(va=va,vb=vb,vab=vab,out=sumof[frozenset((s2,pab))]))
    return res if len(res)==2 else None
ro=outwire(ROOT)
print('root mux coords:',ro)
NODE[0]=ROOT; OUT[0]=ro   # node id 0 == the root combiner
tree[0]=(RA,RB)

# ---- slot links: free wire w  <->  defined z, via atom  M*(w - z) - handle ----
linkre=[re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]'), re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-')]
link={}   # free wire w -> z
for a in E.res:
    m=linkre[0].match(a)
    if m: u,z=int(m.group(1)),int(m.group(2))
    else:
        m=linkre[1].match(a)
        if not m: continue
        u,z=int(m.group(2)),int(m.group(3))
    if u not in defrhs and z in defrhs: link[u]=z
    elif z not in defrhs and u in defrhs: link[z]=u
print('slot links found:',len(link))

# consistency: parent's value wire == link of child's out
good=0; bad=[]
for n,N in NODE.items():
    for side,ch in (('va',N['a']),('vb',N['b'])):
        if tree[ch] is None: continue
        cout={d['out'] for d in OUT[ch]}
        pv={d[side] for d in OUT[n]}
        if {link.get(w) for w in pv}==cout: good+=1
        else: bad.append((n,ch,sorted(pv),sorted(cout),[link.get(w) for w in pv]))
print('parent-slot <-> child-out link consistency: good %d  bad %d'%(good,len(bad)))
for b in bad[:3]: print('   ',b)
pickle.dump(dict(NODE=NODE,OUT=OUT,tree=tree,live=live,dead=dead,link=link),open('model.pkl','wb'))
