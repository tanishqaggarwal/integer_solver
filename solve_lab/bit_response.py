import json
import heal_harness as H
from collections import defaultdict
p=H.p
pinrec=json.load(open('pinrec.json'))
# group pins by selector
bysel=defaultdict(list)
for i,sel,tgt,const,coef,handle in pinrec: bysel[sel].append((tgt,const%p,handle))
bits=sorted(bysel)
print(f"{len(bits)} control bits")
# baseline: fc_partial (new quadrant, x_15298=1)
d=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
LOADS=[11150,25739,37758]  # L1,L2,L3
S,T=35389,6671
base={n:V[n]%p for n in LOADS+[S,T,15298]}
print(f"baseline: x_15298={V[15298]%p}, L=[{[base[n] for n in LOADS]}]")
print(f"  S%p={base[S]}, T%p={base[T]}")
# desc for incremental
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
# single-bit activation: set bit=1, targets=const, measure loads
resp={}
for b in bits:
    saved={}
    saved[b]=V[b]; V[b]=1
    for tgt,const,handle in bysel[b]:
        saved[tgt]=V[tgt]; V[tgt]=const
    # forward affected
    aff=set(desc_of[b])
    for tgt,const,handle in bysel[b]: aff|=set(desc_of[tgt])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
    dL=tuple((V[n]-base_prev)%p if False else (V[n]%p) for n,base_prev in [(x,base[x]) for x in LOADS])
    dS=(V[S]%p); dT=(V[T]%p); dMux=(V[15298]%p)
    resp[b]=(dMux, dS, dT)
    # restore
    for v,x in saved.items(): V[v]=x
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
# analyze: how many bits change S or T?
nz=[(b,r) for b,r in resp.items() if r[1]!=base[S] or r[2]!=base[T]]
muxch=[(b,r[0]) for b,r in resp.items() if r[0]!=base[15298]]
print(f"\nbits that change S or T when activated: {len(nz)}")
print(f"bits that change x_15298: {len(muxch)}: {[b for b,_ in muxch][:20]}")
for b,r in nz[:20]:
    print(f"  bit x_{b}: mux={r[0]}, S={r[1]}, T={r[2]}")
json.dump({str(b):list(r) for b,r in resp.items()}, open('bit_resp.json','w'))
