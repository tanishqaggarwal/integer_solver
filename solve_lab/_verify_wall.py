import heal_harness as H
p=H.p
# gadget_handled: gadget satisfied -> x_2239 should be ≡0 mod p, x_31731=0
d=H.loadd('gadget_handled.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); V=H.val
print("gadget_handled (gadget SAT, G1/G2 open):")
print(f"  x_2239 mod p = {V[2239]%p}  (==0? {V[2239]%p==0})")
print(f"  x_31731 = {V[31731]}  (==0? {V[31731]==0})")
print(f"  x_9106 mod (13523997*p) = {V[9106]%(13523997*p)}  (==0? {V[9106]%(13523997*p)==0})")
print(f"  x_9118 mod p = {V[9118]%p}")
print(f"  x_7068 mod p = {V[7068]%p}   (x_9118==x_7068 mod p? {V[9118]%p==V[7068]%p})")
# my G1=G2=0 config: x_2239 mod p should be nonzero
d2=H.loadd('best_agentA_39022.json')
base={v:d2.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1; base[9118]=base[7068]; base[8731]=base[4432]
for v in H.freeinp: H.val[v]=base[v]
H.forward(); V=H.val
print("\nMy G1=G2=0 config (x_9118=x_7068, x_8731=x_4432):")
print(f"  G1={7376877*V[642]+V[2099]-V[7068]}, G2={V[4432]-V[19964]-V[28730]}")
print(f"  x_2239 mod p = {V[2239]%p}  (==0? {V[2239]%p==0})")
print(f"  -> gadget REQUIRES x_2239≡0 mod p but it's pinned nonzero => WALL confirmed")
# The residue x_2239 mod p as function of x_7068 mod p, x_4432 mod p: can moving x_7068 fix it?
# check: does x_2239 mod p change if I move x_7068 by 1 (keeping x_9118=x_7068)?
from collections import defaultdict
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
t7068=sorted(set(desc_of[7068])|set(desc_of[9118]))
V[7068]=base[7068]+1; V[9118]=base[7068]+1
for k in t7068: V[H.order[k]]=eval(H.gcode[k],ns)
print(f"\n  after x_7068+=1 (x_9118 tracks): x_2239 mod p = {V[2239]%p} (changed? {V[2239]%p!=83965570201604323522827503914887430882251953270432913416137551123691745729744})")
