import json
import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('sy_regime11_39018.json')
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
o8,o9=V[8731],V[9118]
def q27_43():
    return V[27177],V[4306]
b27,b43=q27_43()
# exact slopes (verified linear)
V[9118]=o9+1; fwd_from([9118]); s9_27=V[27177]-b27; s9_43=V[4306]-b43; V[9118]=o9; fwd_from([9118])
V[8731]=o8+1; fwd_from([8731]); s8_27=V[27177]-b27; s8_43=V[4306]-b43; V[8731]=o8; fwd_from([8731])
# verify linearity with step 1000
V[9118]=o9+1000; fwd_from([9118]); chk=V[27177]-b27; V[9118]=o9; fwd_from([9118])
print(f"linear x_27177 in x_9118: {chk==1000*s9_27}")
print(f"s8_27={s8_27} s9_27={s9_27}")
print(f"s8_43={s8_43} s9_43={s9_43}")
# x_27177 = b27 + s8_27*u8 + s9_27*u9  ; x_4306 = b43 + s8_43*u8 + s9_43*u9
# Cond(1): 15964591*x_27177 + 13881285*x_4306 = 0
A=15964591; B=13881285
# = (A*b27+B*b43) + (A*s8_27+B*s8_43)*u8 + (A*s9_27+B*s9_43)*u9 = 0
C1=A*b27+B*b43; P8=A*s8_27+B*s8_43; P9=A*s9_27+B*s9_43
print(f"\nCond1: {P8}*u8 + {P9}*u9 = {-C1}")
import math
g=math.gcd(P8,P9)
print(f"gcd(P8,P9)={g}, divides C1? {(-C1)%g==0}")
