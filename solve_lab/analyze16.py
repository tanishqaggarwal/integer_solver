import json
import heal_harness as H
from collections import defaultdict
p=H.p
# rebuild g1g2_closed state
d=H.loadd('gadget_handled.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
V[7068]=7376877*V[642]+V[2099]; V[4432]=V[19964]+V[28730]
fwd_from([7068,4432])
F=set(H.fails())
print(f"fails after G1G2 close: {len(F)}: {sorted(F)}")
json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('g1g2_closed.json','w'))
# footprint analysis (structural)
eq_frees=[]
for i in range(len(H.lines)):
    s=set()
    for v in H.eqvars[i]: s|=H.anc.get(v, {v} if v in H.freeinp else set())
    eq_frees.append(s & H.freeinp)
foot=defaultdict(set)
for i,fs in enumerate(eq_frees):
    for f in fs: foot[f].add(i)
broken=F
brk_frees=set()
for i in broken: brk_frees|=eq_frees[i]
# classify each free by how many CONSTRAINING (non-currently-satisfied-def) eqs it touches outside broken
# Use forward-based ripple: perturb each candidate free, count new fails
Fset=set(H.fails())
cand=sorted(brk_frees)
print(f"\n{len(cand)} frees feed the 16 broken eqs. Measuring real ripple (forward-based)...")
results=[]
for w in cand:
    old=V[w]; V[w]=old+1; fwd_from([w])
    nf=set(H.fails())
    outside=len(nf-broken-Fset)  # new breaks outside the 16
    helps=len(broken-nf)         # how many of 16 it can affect
    V[w]=old; fwd_from([w])
    results.append((w,outside,helps))
results.sort(key=lambda r:(r[1],-r[2]))
print("free: (new-breaks-outside-16, #of-16-affected)")
for w,o,h in results[:30]:
    print(f"  x_{w}: outside={o} affects16={h} {'PRIVATE' if o==0 else ''}")
