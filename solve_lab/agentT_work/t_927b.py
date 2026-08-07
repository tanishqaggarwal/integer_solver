#!/usr/bin/env python3
"""AUDIT T17b -- close the independence loop on the 927.
t_927.py borrowed L's 3,681-cofactor list.  Derive the handle family from F's parse ALONE:
every atom of shape (h - (P*u)) where P is a wire whose deliverable value is exactly p."""
import os,sys,pickle,collections,re,json
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
v2a=collections.defaultdict(list)
for i,a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
pat=re.compile(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$')
fam=[]
for i,a in enumerate(names):
    m=pat.match(a.replace(' ',''))
    if m: fam.append((i,int(m.group(1)),int(m.group(2)),int(m.group(3))))
print('atoms of shape (h-(P*u)) in F\'s parse: %d'%len(fam))
Pcount=collections.Counter(x[2] for x in fam)
print('distinct P-wires used: %d ; commonest: %s'%(len(Pcount),Pcount.most_common(5)))
H=pickle.load(open(os.path.join(LAB,'agentL_work','handles.pkl'),'rb'))
Lu=set(H['handle'])
myu={x[3] for x in fam}
print('\ncofactor set derived from F alone : %d'%len(myu))
print('agent L\'s cofactor list           : %d'%len(Lu))
print('identical sets? %s   (F-only minus L: %d, L minus F-only: %d)'%(myu==Lu,len(myu-Lu),len(Lu-myu)))
cv=collections.Counter()
for i,h,Pv,u in fam:
    others=[j for j in v2a[h] if j!=i]
    if len(others)!=1: cv['NOGUARD']+=1; continue
    g=names[others[0]].replace(' ','')
    mm=re.search(r'\((\d+)\*x%d(?![0-9])\)'%h,g)
    cv[int(mm.group(1)) if mm else 1]+=1
c1=cv[1]; cg=sum(n for c,n in cv.items() if isinstance(c,int) and c>1)
print('\n== 927 rebuilt with NOTHING borrowed ==')
print('   c == 1 : %d'%c1)
print('   c >  1 : %d   <- the 927'%cg)
print('   distinct c>1 values: %d (each handle its own multiplier? %s)'%(
      len([c for c in cv if isinstance(c,int) and c>1]), len([c for c in cv if isinstance(c,int) and c>1])==cg))
