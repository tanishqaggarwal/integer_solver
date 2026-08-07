#!/usr/bin/env python3
"""AUDIT T21 -- (a) finish the 26 reconciliation, (b) what do the 278 multi-hop aliases route through?"""
import os,sys,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms); avars=[frozenset(vars_of(atoms[a])) for a in names]
v2a=collections.defaultdict(set)
for i,vs in enumerate(avars):
    for u in vs: v2a[u].add(i)
par={}
def find(x):
    par.setdefault(x,x)
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
def uni(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: par[ra]=rb
cp=re.compile(r'^\(x(\d+)-x(\d+)\)$')
for a in names:
    m=cp.match(a.replace(' ',''))
    if m: uni(int(m.group(1)),int(m.group(2)))
PCLASS={x for x in par if find(x)==find(26064)}
print('== (a) the 7 L counts and I missed: is the p-wire the SECOND factor? ==')
for h,u in [(14163,13271),(9145,18624),(1937,23452),(32117,24033),(30351,24874),(14262,31574),(21953,34951)]:
    a=[names[i] for i in v2a[u] if names[i].replace(' ','').startswith('(x%d-'%h)]
    m=re.match(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$',a[0].replace(' ','')) if a else None
    if m:
        f1,f2=int(m.group(2)),int(m.group(3))
        print('   %-34s first=x%-6d inPclass=%-5s   second=x%-6d inPclass=%s'%(a[0][:34],f1,f1 in PCLASS,f2,f2 in PCLASS))
print('   -> if second-factor is the p-wire, these are the SAME family with operands swapped.')
print('\n== (a2) do any of the 33 atoms L omits touch the 25 baseline-failing equations? ==')
import checker as CK
codes,varsets=CK.load_equations()
v2e=collections.defaultdict(set)
for e,vs in enumerate(varsets):
    for u in vs: v2e[u].add(e)
M25=[5324,9041,11226,15558,21000,22534,22997,28929,29330,32026,35512,38051,
     12231,12270,12350,14584,18673,22044,29125]
EXTRA=[242,1997,2922,4671,9701,10422,10462,13632,15035,15120,16206,17073,17286,20997,21170,
       21728,23504,23750,24511,25341,25530,27848,28640,28817,28838,29966,31357,31676,32767,
       33098,34831,35531,36342]
inc=[u for u in EXTRA if v2e[u]&set(M25)]
print('   of the 33 omitted cofactors, incident to the baseline-failing set: %d  %s'%(len(inc),inc))
print('\n== (b) the 278 multi-hop aliases ==')
Mo=pickle.load(open(os.path.join(LAB,'agentL_work','full_model.pkl'),'rb'))
OUT=Mo['OUT']; tree=Mo['tree']
pairs=[]
for n in tree:
    if n not in OUT or not tree[n] or len(tree[n])!=2: continue
    ca,cb=tree[n]
    for j,slot in enumerate(OUT[n]):
        for side,ch in (('va',ca),('vb',cb)):
            w=slot.get(side)
            if ch not in OUT or w is None: continue
            co=OUT[ch][j].get('out')
            if co is None or w==co: continue
            pairs.append((w,co))
multi=[(w,co) for w,co in pairs if not (v2a[w]&v2a[co])]
print('   aliased %d, one-atom %d, multi-hop %d'%(len(pairs),len(pairs)-len(multi),len(multi)))
# shortest atom-path from w to co
def path(w,co,maxd=4):
    seen={w}; fr=[(w,[])]
    for dep in range(maxd):
        nf=[]
        for x,pa in fr:
            for i in v2a[x]:
                for y in avars[i]:
                    if y in seen: continue
                    if y==co: return pa+[i]
                    seen.add(y); nf.append((y,pa+[i]))
        fr=nf
        if not fr: break
    return None
lens=collections.Counter(); shapes=collections.Counter(); viaP=0; nopath=0
for w,co in multi:
    pth=path(w,co)
    if pth is None: nopath+=1; continue
    lens[len(pth)]+=1
    hit=any(avars[i]&PCLASS for i in pth)
    if hit: viaP+=1
    shapes[' | '.join(re.sub(r'x\d+','x',names[i].replace(' ','')) for i in pth)]+=1
print('   shortest-path length histogram: %s   (no path <=4: %d)'%(dict(lens),nopath))
print('   multi-hop pairs whose PATH touches a p-class wire: %d of %d'%(viaP,len(multi)-nopath))
print('   commonest path shapes:')
for s,n in shapes.most_common(8): print('      %-74s %d'%(s[:74],n))
