#!/usr/bin/env python3
"""AUDIT T23 -- CLOSE loose end 1.  The 278 'multi-hop' aliases look like an artifact of my own
index-matched pairing (OUT[n][j] <-> OUT[child][j]).  L's calib2 recovered a per-node parent/child
COORDINATE ALIGNMENT (188 orient=1, 67 orient=0), so for flipped nodes the correct partner is
OUT[child][1-j].  Re-pair allowing the cross and re-measure."""
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
for a in names:
    m=re.match(r'^\(x(\d+)-x(\d+)\)$',a.replace(' ',''))
    if m: uni(int(m.group(1)),int(m.group(2)))
PCLASS={x for x in par if find(x)==find(26064)}
prod=re.compile(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$')
HANDLE={}      # h -> u   for h = P*u with P in p-class (either operand order)
for i,a in enumerate(names):
    m=prod.match(a.replace(' ',''))
    if not m: continue
    h,f1,f2=int(m.group(1)),int(m.group(2)),int(m.group(3))
    if f1 in PCLASS: HANDLE[h]=f2
    elif f2 in PCLASS: HANDLE[h]=f1
print('handles h = p*u (either operand order): %d'%len(HANDLE))
Mo=pickle.load(open(os.path.join(LAB,'agentL_work','full_model.pkl'),'rb'))
OUT=Mo['OUT']; tree=Mo['tree']
same=cross=neither=0; slackP=0; slackOther=[]
for n in tree:
    if n not in OUT or not tree[n] or len(tree[n])!=2: continue
    ca,cb=tree[n]
    for j,slot in enumerate(OUT[n]):
        for side,ch in (('va',ca),('vb',cb)):
            w=slot.get(side)
            if ch not in OUT or w is None: continue
            o_same=OUT[ch][j].get('out'); o_cross=OUT[ch][1-j].get('out')
            hit=None
            if o_same is not None and (w==o_same or (v2a[w]&v2a[o_same])): same+=1; hit=o_same
            elif o_cross is not None and (w==o_cross or (v2a[w]&v2a[o_cross])): cross+=1; hit=o_cross
            else: neither+=1; continue
            # the alias atom's third variable = the slack; is it a p-handle?
            for i in v2a[w]&v2a[hit]:
                extra=avars[i]-{w,hit}
                for s in extra:
                    if s in HANDLE: slackP+=1
                    else: slackOther.append((n,side,s))
                break
tot=same+cross+neither
print('\n== parent/child links re-paired allowing the coordinate cross ==')
print('   aliased via SAME coordinate index  : %d'%same)
print('   aliased via CROSSED index          : %d'%cross)
print('   still no one-atom alias            : %d'%neither)
print('   TOTAL                              : %d'%tot)
print('\n== the slack wire in each alias atom ==')
print('   slack is a p-handle (= p*u)  : %d'%slackP)
print('   slack is something else      : %d  %s'%(len(slackOther),slackOther[:4]))
