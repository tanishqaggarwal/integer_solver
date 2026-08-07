"""FULL model over the TRUE global OR tree (root x9274, 384 leaves: 256 free + 128 dead)."""
import sys, os, collections, pickle, re, time, json
F='/home/user/integer_solver/solve_lab/agentF_work'; sys.path.insert(0,F)
from fwd import Engine, NV
from parse import node_str
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
H=pickle.load(open('handles.pkl','rb')); appearP=H['appearP']
ORS=pickle.load(open('ors.pkl','rb')); ROOT=9274
resby=collections.defaultdict(list)
for a in E.res:
    for u in vars_of(E.atoms[a]): resby[u].append(a)
def deref(v):
    seen=set()
    while v in defrhs and defrhs[v][0]=='v' and v not in seen: seen.add(v); v=defrhs[v][1]
    return v
prodof=collections.defaultdict(list); notof=collections.defaultdict(list)
gated=collections.defaultdict(list); sumof={}
for w,r in defrhs.items():
    if r[0]=='*' and r[1][0]=='v' and r[2][0]=='v':
        a,b=deref(r[1][1]),deref(r[2][1]); prodof[frozenset((a,b))].append(deref(w))
        gated[a].append((deref(w),deref(r[2][1]) if deref(r[1][1])==a else deref(r[1][1])))
        gated[b].append((deref(w),deref(r[1][1]) if deref(r[2][1])==b else deref(r[2][1])))
    elif r[0]=='-' and r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': notof[deref(r[2][1])].append(deref(w))
    elif r[0]=='+' and r[1][0]=='v' and r[2][0]=='v': sumof[frozenset((deref(r[1][1]),deref(r[2][1])))]=deref(w)
NODE={}; OUT={}
for n,(a,b) in ORS.items():
    sa=[x for nb in notof[b] for x in prodof.get(frozenset((a,nb)),[])]
    sb=[x for na in notof[a] for x in prodof.get(frozenset((b,na)),[])]
    sab=[s for s in prodof.get(frozenset((a,b)),[]) if any(w!=n for w,_ in gated.get(s,[]))]
    N=dict(a=a,b=b,sa=sa,sb=sb,sab=sab,
           ga=[g for s in sa for g in gated[s]],gb=[g for s in sb for g in gated[s]],
           gab=[g for s in sab for g in gated[s]])
    res=[]
    for (pa,va) in N['ga']:
        c=[(pb,vb) for (pb,vb) in N['gb'] if frozenset((pa,pb)) in sumof]
        if len(c)!=1: res=None; break
        pb,vb=c[0]; s2=sumof[frozenset((pa,pb))]
        c2=[(pab,vab) for (pab,vab) in N['gab'] if frozenset((s2,pab)) in sumof]
        if len(c2)!=1: res=None; break
        pab,vab=c2[0]; res.append(dict(va=va,vb=vb,vab=vab,out=sumof[frozenset((s2,pab))]))
    NODE[n]=N; OUT[n]=res if (res and len(res)==2) else None
print('nodes %d ; with clean mux %d ; without %s'%(len(NODE),sum(1 for x in OUT.values() if x),
      [n for n,x in OUT.items() if not x][:6]))
tree={n:ORS[n] for n in ORS}
leaves=set()
def col(n):
    if n not in ORS: leaves.add(n); return
    col(ORS[n][0]); col(ORS[n][1])
col(ROOT)
for L in leaves: tree[L]=None
live=[L for L in leaves if L not in defrhs]; dead=[L for L in leaves if L in defrhs]
print('leaves %d live %d dead %d ; all dead literal-0: %s'%(len(leaves),len(live),len(dead),
      all(defrhs[d]==('c',0) for d in dead)))
# links
link={}
lr=[re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]'), re.compile(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-')]
for a in E.res:
    m=lr[0].match(a)
    if m: u,z=int(m.group(1)),int(m.group(2))
    else:
        m=lr[1].match(a)
        if not m: continue
        u,z=int(m.group(2)),int(m.group(3))
    if u not in defrhs and z in defrhs: link[u]=z
    elif z not in defrhs and u in defrhs: link[z]=u
# subtree leaves
sub={}
def subl(n):
    if n in sub: return sub[n]
    sub[n]=[n] if tree[n] is None else subl(tree[n][0])+subl(tree[n][1])
    return sub[n]
subl(ROOT)
order=[]
def post(n):
    if tree[n] is None: return
    for c in tree[n]: post(c)
    order.append(n)
post(ROOT)
def run(v):
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    return vv,E.run(vv)
# ---- leaf pins (numeric) ----
leafnode={}
for n in NODE:
    for side,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])):
        if tree[ch] is None: leafnode[ch]=(n,side)
t0=time.time(); PIN={}
for L in live:
    n,side=leafnode[L]; ws=[d[side] for d in OUT[n]]
    _,r0=run({L:1}); Cs=[]
    for w in ws:
        _,r1=run({L:1,w:1}); cand=None
        for a in appearP.get(w,[]):
            i=E.residx[a]; f0=r0[i]%p; sl=(r1[i]-r0[i])%p
            if f0 and sl:
                c=(-f0)*pow(sl,p-2,p)%p
                cand=c if cand is None else ('CONFLICT' if cand!=c else cand)
        Cs.append(cand)
    PIN[L]=(ws,Cs)
badpin=[L for L,(w,c) in PIN.items() if any(x is None or x=='CONFLICT' for x in c)]
print('pins: %d leaves, bad %d %s (%.0fs)'%(len(PIN),len(badpin),badpin[:5],time.time()-t0))
pickle.dump(dict(NODE=NODE,OUT=OUT,tree=tree,live=live,dead=dead,link=link,sub=sub,order=order,
                 PIN=PIN,ROOT=ROOT,leafnode=leafnode),open('full_model.pkl','wb'))
