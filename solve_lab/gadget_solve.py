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
print(f"baseline fails: {len(F0)}: {sorted(F0)}")
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
def getq(): return V[27177],V[4306]
# check linearity of (x_27177,x_4306) in x_8731,x_9118
o8,o9=V[8731],V[9118]
base=getq()
V[8731]=o8+1; fwd_from([8731]); d8=getq(); V[8731]=o8; fwd_from([8731])
V[9118]=o9+1; fwd_from([9118]); d9=getq(); V[9118]=o9; fwd_from([9118])
V[8731]=o8+2; fwd_from([8731]); d8b=getq(); V[8731]=o8; fwd_from([8731])
# slopes
a=(d8[0]-base[0]); b=(d9[0]-base[0])  # d x_27177
c=(d8[1]-base[1]); e=(d9[1]-base[1])  # d x_4306
lin8 = (d8b[0]-base[0])==2*a  # x_27177 linear in x_8731?
print(f"linearity check (x_27177 vs x_8731): {lin8}")
print(f"slopes: dx27177/d8731={a}, dx27177/d9118={b}")
print(f"        dx4306/d8731={c}, dx4306/d9118={e}")
print(f"base x_27177={base[0]}\n     x_4306={base[1]}")
# solve linear system: base + a*du8 + b*du9 = 0 ; base1 + c*du8 + e*du9 = 0
det=a*e-b*c
print(f"det={det}")
if det!=0:
    # Cramer over rationals -> need integer. Solve exactly.
    from fractions import Fraction
    du8=Fraction(-base[0]*e + b*base[1], det)
    du9=Fraction(-a*base[1] + c*base[0], det)
    print(f"du8={du8} (int={du8.denominator==1}), du9={du9} (int={du9.denominator==1})")
    if du8.denominator==1 and du9.denominator==1:
        V[8731]=o8+int(du8); V[9118]=o9+int(du9)
        fwd_from([8731,9118])
        print(f"after solve: x_27177={V[27177]}, x_4306={V[4306]}")
        F=set(H.fails())
        print(f"fails: {len(F)}  fixed:{sorted(F0-F)}  broken:{sorted(F-F0)}")
