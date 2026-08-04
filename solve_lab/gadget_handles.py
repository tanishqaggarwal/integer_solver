import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('gadget_zeroed.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
# x_6947 = 6122989 * (x_2239//p)
M=V[2239]//p
V[6947]=6122989*M
# x_33168=0 already
V[33168]=0
# x_9629 = x_9106/13523997 via x_950 (x_9629 = x_30095*x_950)
Q=V[9106]//13523997
x30095=V[30095]
print(f"x_30095={x30095}, Q%x_30095=={Q%x30095 if x30095 else 'div0'}")
if x30095!=0 and Q%x30095==0:
    V[950]=Q//x30095
else:
    print("  x_30095 does not divide Q; need alt handle")
fwd_from([6947,33168,950])
print(f"atom17897 = x_9106-13523997*x_9629 = {V[9106]-13523997*V[9629]}")
print(f"atom20866 = 6122989*x_2239 - x_23754 = {6122989*V[2239]-V[23754]}")
print(f"atom20868 = x_21279*x_31731 + x_35619 = {V[21279]*V[31731]+V[35619]}")
print(f"atom34232 = x_21279*x_31731 = {V[21279]*V[31731]}")
F=set(H.fails())
print(f"\nFAILS: {len(F)}: {sorted(F)}")
print(f"fixed: {sorted(F0-F)}")
print(f"broken: {sorted(F-F0)}")
json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('gadget_handled.json','w'))
