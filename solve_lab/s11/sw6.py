"""Full construction with the x4287 bit ON.

x4287 = 1  =>  x21279 = 1  =>  x7075 = 0, so the checkpoint's two hard congruences
(p | x9118, p | x8731) disappear.  In exchange the x21279 channel switches on three new ones:

    a22233 : 6122989*x2239*x21279 = x23754 = p*x6947    =>  p | x2239
    a22235 : x21279*x31731 = -x35619 = -p*x33168        =>  p | x31731
    a19088 : x9106*x21279 = 13523997*x9629 = ...*p      =>  p | x9106

and    x2239   = 3494591*x27177 + 14240157*x4306
       x31731  = 15964591*x27177 + 13881285*x4306
so all three reduce to  x27177 == 0  and  x4306 == 0  (mod p), and

    x27177 = x17925^2*(x9118 + x31861 + x6418 + x24453) - x27019^2      -- affine in x9118
    x4306  = (x8731 + x14865)*x17925 - x27019*(x31861 - x9118)          -- affine in x8731

Two knobs, two congruences.  Solve them mod p by fitting each response at two points.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; Q=7376877
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def rep(v,tag,show=20):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken={len(B)} {B[:show]} failing={len(F)} score={L.NEQ-len(F)}")
    return B,F
BIG1=[abs(c) for m,c in L.polys[3568].items() if m==(4287,)][0]
BIG2=[abs(c) for m,c in L.polys[3570].items() if m==(4287,)][0]
x7068_0=v[7068]
L.ripple(v,{4287:1, 31861:BIG1, 14865:BIG2})
print(f"bit on: x7075={v[7075]} x21279={v[21279]}")

def probe(v, knob, target, t0, t1):
    a=list(v); L.ripple(a,{knob:t0}); y0=a[target]%P
    b=list(v); L.ripple(b,{knob:t1}); y1=b[target]%P
    B=(y1-y0)*pow((t1-t0)%P,-1,P)%P
    return y0, B, t0
# 1) x9118 kills x27177
y0,Bc,t0 = probe(v, 9118, 27177, 0, 1)
assert Bc, "x27177 not affine-responsive to x9118"
w9 = (t0 - y0*pow(Bc,-1,P)) % P
# integer lift of x9118 close to x7068_0 (keeps x7068 = x2099 + Q*x642 near its old value)
n9 = w9 + P*((x7068_0-w9)//P)
# 2) x8731 kills x4306, with x9118 already at n9
v1=list(v); L.ripple(v1,{9118:n9})
y0,Bd,t0 = probe(v1, 8731, 4306, 0, 1)
assert Bd, "x4306 not affine-responsive to x8731"
n8 = (t0 - y0*pow(Bd,-1,P)) % P
seeds={9118:n9, 8731:n8, 29854:0, 1329:0, 31864:0, 10903:0,
       17325:0, 642:0, 9413:v[28730]//P, 28730:P*(v[28730]//P), 21574:0, 1844:0}
L.ripple(v,seeds)
print(f"x27177 %p == 0 : {v[27177]%P==0}    x4306 %p == 0 : {v[4306]%P==0}")
print(f"x2239 %p==0: {v[2239]%P==0}   x31731 %p==0: {v[31731]%P==0}   x9106 %p==0: {v[9106]%P==0}")
# now set the three quotient handles the new checks want
h={}
if v[2239]%P==0:  h[6947]  = 6122989*v[2239]//P
if v[31731]%P==0: h[33168] = -v[31731]//P
if v[9106]%P==0:  h[950]   = v[9106]//(13523997*P) if v[9106]%(13523997*P)==0 else None
print("handle for x950 exact:", h.get(950) is not None, " x9106 %(13523997*p)==0:", v[9106]%(13523997*P)==0)
h={k:x for k,x in h.items() if x is not None}
L.ripple(v,h)
B,F=rep(v,'after full construction')
for a in B[:20]: print(f"   a{a} out={str(L.atom_out.get(a)):12s} in {len(L.atom2eq.get(a,{}))} eqs")
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','sw6_out.json'),'w'))
