#!/usr/bin/env python3
"""Agent P: wiring graph of the 383 law-blocks."""
import pickle,sys,json
from collections import Counter,defaultdict,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
B=pickle.load(open(W+'blocks.pkl','rb'))
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))
print("blocks:",len(B))

# each block: input pairs X=(i2,i3) Y=(i1,i4) ; law output Z=(i5,i6) ; mux out
# mux terms: list of (out,k,u,v) triples for each of the two coordinate muxes
def muxsrc(b):
    """for each of 2 coord muxes, the 3 source vars in order of the 3 gate products"""
    r=[]
    for terms in b['mux']:
        srcs=[]
        for (o,k,u,v) in terms:
            srcs.append((u,v))
        r.append(srcs)
    return r

# identify liveness vars a,b
live=[b['live'] for b in B]
# build var -> producing block (mux outputs)
prod_of={}
for i,b in enumerate(B):
    for o in b['muxout']: prod_of[o]=i
print("mux output vars:",len(prod_of))

# leaf coordinate vars
leafvar={}
for si,(a,bb,k) in enumerate(LEAVES): leafvar.setdefault(a,[]).append(bb)

edges=defaultdict(set); srckind=Counter()
inputs_of=[]
for i,b in enumerate(B):
    xs=[b['i1'],b['i2'],b['i3'],b['i4']]
    src=[]
    for x in xs:
        if x in prod_of: src.append(('B',prod_of[x])); edges[prod_of[x]].add(i)
        else: src.append(('?',x))
    inputs_of.append(src)
    srckind[tuple(k for k,_ in src)]+=1
print("input source kinds:",srckind.most_common(6))

allleafcoords=set(x for v in leafvar.values() for x in v)
kinds=Counter()
for i,b in enumerate(B):
    for x in (b['i1'],b['i2'],b['i3'],b['i4']):
        kinds['block' if x in prod_of else ('leafcoord' if x in allleafcoords else 'other')]+=1
print("input var kinds:",kinds)
pickle.dump({'prod_of':prod_of,'inputs_of':inputs_of,'edges':{k:sorted(v) for k,v in edges.items()}},open(W+'wire.pkl','wb'))
