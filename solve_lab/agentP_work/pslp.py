#!/usr/bin/env python3
"""Agent P: build & validate the SLP forward evaluator."""
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb'))
rows,AP=D['rows'],D['AP']
O=pickle.load(open(W+'order.pkl','rb'))
topo=O['topo']; C=O['C']
NV=38748
pos={a:i for i,a in enumerate(topo)}

# first appearance
first={}
for i,a in enumerate(topo):
    for m in AP[a]:
        for x in m:
            if x not in first: first[x]=i
byfirst=defaultdict(list)
for x,i in first.items(): byfirst[i].append(x)
print("vars first-appearing per atom:",sorted(Counter(len(v) for v in byfirst.values()).items()))

outof=[-1]*len(AP); definedat={}
amb=0; noout=0
for i,a in enumerate(topo):
    new=[x for x in byfirst.get(i,[]) if x in C[a]]
    if len(new)==1:
        outof[a]=new[0]; definedat[new[0]]=i
    elif len(new)==0:
        noout+=1
    else:
        amb+=1
print("atoms with unique new-candidate output:",sum(1 for o in outof if o>=0),"no-output:",noout,"ambiguous:",amb)
free=[x for x in range(NV) if x not in definedat]
print("free vars:",len(free))

# forward evaluate
v=[None]*NV
gold=[0]*NV
for k,val in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items():
    gold[int(k[2:])]=int(val)
for x in free: v[x]=gold[x]

bad=[];unres=[]
for i,a in enumerate(topo):
    o=outof[a]
    if o<0: continue
    ap=AP[a]
    # solve  co*o + rest = 0
    co=0; rest=0; ok=True
    for m,c in ap.items():
        if m==(o,): co+=c; continue
        t=c
        for x in m:
            if v[x] is None: ok=False; break
            t*=v[x]
        if not ok: break
        rest+=t
    if not ok: unres.append(a); continue
    if co==0 or rest % co: bad.append((a,co)); v[o]=None; continue
    v[o]=-rest//co
nones=sum(1 for x in v if x is None)
print("unresolvable (missing input):",len(unres),"non-divisible:",len(bad),"unset:",nones)
mism=[x for x in range(NV) if v[x] is not None and v[x]!=gold[x]]
print("mismatch vs deliverable:",len(mism), mism[:10])
pickle.dump({'outof':outof,'topo':topo,'free':free,'definedat':definedat},open(W+'slp.pkl','wb'))
