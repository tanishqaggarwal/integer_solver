#!/usr/bin/env python3
"""Agent P: wiring via mod-P equality classes."""
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
g=[0]*NV
for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): g[int(k[2:])]=int(v)
PV={x for x in range(NV) if g[x]==P}

par=list(range(NV))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def uni(a,b):
    a,b=find(a),find(b)
    if a!=b: par[a]=b

neq=0
for ap in AP:
    d=defaultdict(int)
    for m,c in ap.items():
        if any(x in PV for x in m): continue
        d[m]+=c
    d={m:c for m,c in d.items() if c%P}
    if len(d)==2:
        (m1,c1),(m2,c2)=list(d.items())
        if len(m1)==1 and len(m2)==1 and (c1+c2)%P==0:
            uni(m1[0],m2[0]); neq+=1
print("mod-P equality merges:",neq)
cls=defaultdict(list)
for x in range(NV): cls[find(x)].append(x)
print("classes:",len(cls))
R=lambda x: find(x)

prod_of={}
for i,b in enumerate(B):
    for o in b['muxout']: prod_of[R(o)]=i
print("distinct mux-output classes:",len(prod_of))

leafcoord={}
for si,(a,bb,k) in enumerate(LEAVES): leafcoord[R(bb)]=(a,k)
print("leaf coord classes:",len(leafcoord))

kinds=Counter(); edges=defaultdict(set); srcs=[]
for i,b in enumerate(B):
    ss=[]
    for x in (b['i1'],b['i2'],b['i3'],b['i4']):
        r=R(x)
        if r in prod_of: kinds['block']+=1; ss.append(('B',prod_of[r])); edges[prod_of[r]].add(i)
        elif r in leafcoord: kinds['leaf']+=1; ss.append(('L',leafcoord[r][0]))
        else: kinds['other']+=1; ss.append(('?',x))
    srcs.append(ss)
print("block input kinds:",kinds)
# also the law-output pair (i5,i6): is it the mux 3rd source?
m3=Counter()
for i,b in enumerate(B):
    third=[t[2] for t in b['mux'][0]]+[t[3] for t in b['mux'][0]]
    m3[R(b['i5']) in {R(x) for x in third}]+=1
print("i5 among mux0 sources:",m3)
pickle.dump({'par':[find(x) for x in range(NV)],'prod_of':prod_of,'srcs':srcs,'leafcoord':leafcoord,'edges':{k:sorted(v) for k,v in edges.items()}},open(W+'wire2.pkl','wb'))
