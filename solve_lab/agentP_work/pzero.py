#!/usr/bin/env python3
"""Agent P: variables provably ==0 mod P (handles), then mod-P equality classes + wiring."""
import pickle,sys,json
from collections import Counter,defaultdict,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))
NV=38748
# P-aliases: alias closure from the pinned var
par0=list(range(NV))
def f0(a):
    while par0[a]!=a: par0[a]=par0[par0[a]]; a=par0[a]
    return a
pinP=None
for ap in AP:
    it=list(ap.items())
    if len(it)==2:
        (m1,c1),(m2,c2)=sorted(it,key=lambda z:len(z[0]))
        if len(m1)==1 and len(m2)==1 and c1==-c2:
            a,b=f0(m1[0]),f0(m2[0])
            if a!=b: par0[a]=b
        elif m1==() and len(m2)==1 and abs(c1)==P: pinP=m2[0]
print("P pinned at var:",pinP)
PV={x for x in range(NV) if f0(x)==f0(pinP)}
print("P-alias class size:",len(PV))

Z=set(PV)
var2at=defaultdict(list)
for i,ap in enumerate(AP):
    s=set()
    for m in ap: s.update(m)
    for x in s: var2at[x].append(i)
q=deque(range(len(AP)))
while q:
    i=q.popleft(); ap=AP[i]
    d=defaultdict(int)
    for m,c in ap.items():
        if any(x in Z for x in m): continue
        d[m]=(d[m]+c)%P
    d={m:c for m,c in d.items() if c}
    if len(d)==1:
        m=next(iter(d))
        if len(m)==1 and m[0] not in Z:
            Z.add(m[0])
            for j in var2at[m[0]]: q.append(j)
print("vars ==0 mod P (handles+aliases):",len(Z))

par=list(range(NV))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
neq=0
for ap in AP:
    d=defaultdict(int)
    for m,c in ap.items():
        if any(x in Z for x in m): continue
        d[m]=(d[m]+c)%P
    d={m:c for m,c in d.items() if c}
    if len(d)==2:
        (m1,c1),(m2,c2)=list(d.items())
        if len(m1)==1 and len(m2)==1 and (c1+c2)%P==0:
            a,b=find(m1[0]),find(m2[0])
            if a!=b: par[a]=b; neq+=1
print("mod-P equality merges:",neq)
R=lambda x: find(x)
prod_of={}
for i,b in enumerate(B):
    for j,o in enumerate(b['muxout']): prod_of[R(o)]=(i,j)
leafcoord={}
for (a,bb,k) in LEAVES: leafcoord[R(bb)]=(a,k)
kinds=Counter(); edges=defaultdict(set); srcs=[]
for i,b in enumerate(B):
    ss=[]
    for nm in ('i1','i2','i3','i4','i5','i6'):
        r=R(b[nm])
        if r in prod_of: kinds[nm+'-block']+=1; ss.append(('B',)+prod_of[r]); edges[prod_of[r][0]].add(i)
        elif r in leafcoord: kinds[nm+'-leaf']+=1; ss.append(('L',leafcoord[r][0]))
        else: kinds[nm+'-other']+=1; ss.append(('?',b[nm]))
    srcs.append(ss)
print("block input kinds:",dict(sorted(kinds.items())))
pickle.dump({'par':[find(x) for x in range(NV)],'Z':sorted(Z),'prod_of':prod_of,'srcs':srcs,'leafcoord':leafcoord,
             'edges':{k:sorted(v) for k,v in edges.items()}},open(W+'wire3.pkl','wb'))
