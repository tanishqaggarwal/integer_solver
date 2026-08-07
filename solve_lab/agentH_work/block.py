"""Closure of the free-input <-> check-atom bipartite graph from the residual seed."""
import ev,pickle,json,os
from collections import defaultdict
S=pickle.load(open('support.pkl','rb')); csup=S['csup']; fidx=S['fidx']; inv={i:v for v,i in fidx.items()}
def bits(x):
    o=[]
    while x: b=x&-x; o.append(b.bit_length()-1); x^=b
    return o
# index: free input bit -> checks containing it
f2c=defaultdict(list)
for a,s in csup.items():
    for b in bits(s): f2c[b].append(a)
eq_of=defaultdict(list)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: eq_of[a].append(i)

seed=[688,1618,23001,39067,40608]
F=set(); A=set()
front=set()
for a in seed:
    A.add(a); front.update(bits(csup[a]))
while front:
    nf=set()
    for b in front:
        if b in F: continue
        F.add(b)
        for a in f2c[b]:
            if a not in A:
                A.add(a); nf.update(bits(csup[a]))
    front=nf-F
eqs=set()
for a in A: eqs.update(eq_of[a])
print('CLOSURE: free inputs %d, check atoms %d, equations %d'%(len(F),len(A),len(eqs)))
freevars=sorted(inv[b] for b in F)
json.dump({'free':freevars,'atoms':sorted(A),'eqs':sorted(eqs)},open('block.json','w'))
# how many of the closure's free inputs are boolean-ish?
print('sample free vars:',freevars[:40])
