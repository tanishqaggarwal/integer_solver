import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('gadget_handled.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
print(f"start: {len(F0)} fails: {sorted(F0)}")
# check status of x_7068,x_4432,x_2099,x_19964,x_9118,x_8731
for n in [7068,4432,2099,19964,9118,8731,642,28730,17325,9413]:
    print(f"  x_{n}: {'FREE' if n in H.freeinp else 'gate'} = {str(V[n])[:40]}")
# G1 = 7376877*x_642 + x_2099 - x_7068 ; G2 = x_4432 - x_19964 - x_28730
G1=7376877*V[642]+V[2099]-V[7068]
G2=V[4432]-V[19964]-V[28730]
print(f"G1={G1}\nG2={G2}")
# close: x_7068 free -> x_7068 = 7376877*x_642 + x_2099 ; x_4432 free -> x_4432 = x_19964 + x_28730
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
V[7068]=7376877*V[642]+V[2099]
V[4432]=V[19964]+V[28730]
fwd_from([7068,4432])
G1=7376877*V[642]+V[2099]-V[7068]; G2=V[4432]-V[19964]-V[28730]
print(f"after close: G1={G1}, G2={G2}")
F=set(H.fails())
print(f"FAILS: {len(F)}: {sorted(F)}")
print(f"fixed: {sorted(F0-F)}")
print(f"broken: {sorted(F-F0)}")
if len(F)<11:
    json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('g1g2_closed.json','w'))
    print("saved g1g2_closed.json")
