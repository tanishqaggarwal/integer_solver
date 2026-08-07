"""Transitive free-input support of every variable / check atom / equation, as bitsets."""
import model, pickle, os, time
from collections import Counter, defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
d=model.get(); atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
F=pickle.load(open(os.path.join(HERE,'fwd2.pkl'),'rb'))
tgt=F['tgt']; definer=F['definer']; order=F['order']; free0=F['free0']; checks=F['checks']
NV=38748
fidx={v:i for i,v in enumerate(free0)}
sup=[0]*NV
for v in free0: sup[v]=1<<fidx[v]
t0=time.time()
for v in order:
    a=definer[v]
    s=0
    for u in atom_vars[a]:
        if u!=v: s|=sup[u]
    sup[v]=s
print('support built %.1fs'%(time.time()-t0))
sizes=[bin(sup[v]).count('1') for v in range(NV)]
print('var support size: max',max(sizes),'mean %.1f'%(sum(sizes)/NV))
print('hist(log2 buckets):',Counter(0 if s==0 else s.bit_length() for s in sizes).most_common(20))

# check atom supports
csup={}
for a in checks:
    s=0
    for u in atom_vars[a]: s|=sup[u]
    csup[a]=s
cs=[bin(v).count('1') for v in csup.values()]
print('check support size: max',max(cs),'mean %.1f'%(sum(cs)/len(cs)))
print('check hist:',Counter(min(v,50) for v in cs).most_common(20))
# equation supports (union over its check atoms only? no: over all atoms, but gate atoms are 0 by construction)
esup=[]
for i,(m,sq,tl) in enumerate(eq_terms):
    s=0
    for c,a in tl:
        if a in csup: s|=csup[a]
    esup.append(s)
es=[bin(v).count('1') for v in esup]
print('eq support: max',max(es),'mean %.1f'%(sum(es)/len(es)))
print('eqs with empty check-support:', sum(1 for v in es if v==0))
pickle.dump({'sup':sup,'csup':csup,'esup':esup,'fidx':fidx},open(os.path.join(HERE,'support.pkl'),'wb'))
