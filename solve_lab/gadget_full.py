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
# verify linear-combo relations
r31=V[31731]-(15964591*V[27177]+13881285*V[4306])
r91=V[9106]-(7204959*V[27177]+6822253*V[4306])
r22=V[2239]-(3494591*V[27177]+14240157*V[4306])
print(f"relation residuals: x_31731 {r31}, x_9106 {r91}, x_2239 {r22}")
o8,o9=V[8731],V[9118]
b27,b43=V[27177],V[4306]
V[9118]=o9+1; fwd_from([9118]); s9_27=V[27177]-b27; s9_43=V[4306]-b43; V[9118]=o9; fwd_from([9118])
V[8731]=o8+1; fwd_from([8731]); s8_43=V[4306]-b43; V[8731]=o8; fwd_from([8731])
s8_27=0
# Cond1 exact: A*x27+B*x43=0
A,B=15964591,13881285
C1=A*b27+B*b43; P8=A*s8_43+A*s8_27*B/A if False else A*s8_27+B*s8_43; P9=A*s9_27+B*s9_43
# extended euclid for P8*u8+P9*u9 = -C1
def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)
g,xx,yy=egcd(P8,P9)
assert (-C1)%g==0
u8_0=xx*(-C1//g); u9_0=yy*(-C1//g)
# general: u8=u8_0+ (P9//g)*t ; u9=u9_0-(P8//g)*t
du8=P9//g; du9=-P8//g
# x_27177(t)=b27+s9_27*u9 ; x_4306(t)=b43+s8_43*u8+s9_43*u9
def x27(t): u9=u9_0+du9*t; return b27+s9_27*u9
def x43(t): u8=u8_0+du8*t; u9=u9_0+du9*t; return b43+s8_43*u8+s9_43*u9
# Cond2: 13523997 | 7204959*x27+6822253*x43  ; Cond3: 3494591*x27+14240157*x43 ≡0 mod p
# both are affine in t: compute coefficients
c2_0=(7204959*x27(0)+6822253*x43(0))
c2_1=(7204959*x27(1)+6822253*x43(1))-c2_0
c3_0=(3494591*x27(0)+14240157*x43(0))
c3_1=(3494591*x27(1)+14240157*x43(1))-c3_0
print(f"Cond2: {c2_1}*t + {c2_0} ≡0 mod 13523997")
print(f"Cond3: coef*t+const ≡0 mod p")
# solve Cond2 mod m2, Cond3 mod p, CRT
m2=13523997
# t ≡ t2 mod (m2/gcd)
from math import gcd
def solve_mod(c1,c0,m):
    c1%=m; c0%=m; gg=gcd(c1,m)
    if (-c0)%gg!=0: return None
    m2=m//gg; c1r=(c1//gg)%m2; c0r=((-c0)//gg)%m2
    inv=pow(c1r,-1,m2); return ( (c0r*inv)%m2, m2)
r2=solve_mod(c2_1,c2_0,m2)
r3=solve_mod(c3_1%p,c3_0%p,p)
print(f"t mod branch2: {r2 if r2 else 'NO'}")
print(f"t mod branch3: {(str(r3[0])[:30]+'...',r3[1]) if r3 else 'NO'}")
if r2 and r3:
    # CRT t ≡ r2 mod m2', t ≡ r3 mod p'
    a1,n1=r2; a2,n2=r3
    gg=gcd(n1,n2)
    if (a2-a1)%gg==0:
        lcm=n1//gg*n2
        inv=pow(n1//gg, -1, n2//gg)
        t=(a1 + (n1)*(((a2-a1)//gg*inv)%(n2//gg)))%lcm
        print(f"CRT solution t found (mod {str(lcm)[:20]}...). t={str(t)[:40]}...")
        # apply
        u8=u8_0+du8*t; u9=u9_0+du9*t
        V[8731]=o8+u8; V[9118]=o9+u9
        fwd_from([8731,9118])
        print(f"applied: x_31731={V[31731]}")
        print(f"  x_9106 % 13523997 = {V[9106]%13523997}")
        print(f"  x_2239 % p = {V[2239]%p}")
        F=set(H.fails())
        print(f"FAILS: {len(F)}  fixed:{sorted(F0-F)}\n  broken:{sorted(F-F0)}")
        json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('gadget_zeroed.json','w'))
    else: print("CRT no solution")
