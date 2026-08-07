"""Flip the message bit x4287 on.

a33881: x21279 = x9062*x20434 = x4287*x2081, and a36085: x7075 = 1 - x21279.
So x4287 = 1  =>  x7075 = 0, and BOTH hard congruences evaporate:

    T = 5113045*x7075*x9118 = 0      U = x7075*x8731 = 0

Also x31033 = x20434*x31822 = 0, so x22542 = 0 and x2099 collapses to x25297 = x9118*x21279
= x9118 -- a FREE variable.  Then x7068 = x9118 + 7376877*x642 can be held at its old value
while x642 becomes a multiple of p, so the copy network never moves.

Cost: x4287 is a pinned message bit -- turning it on pins x31861 and x14865.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def rep(v,tag,show=14):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken atoms={len(B)} {B[:show]}  failing={len(F)} score={L.NEQ-len(F)}")
    return B,F
rep(v,'checkpoint')
Q=7376877
BIG1=[abs(c) for m,c in L.polys[3568].items() if m==(4287,)][0]
BIG2=[abs(c) for m,c in L.polys[3570].items() if m==(4287,)][0]
x7068_0=v[7068]
h=int(sys.argv[1]) if len(sys.argv)>1 else 0
seeds={4287:1, 31861:BIG1, 14865:BIG2,
       29854:0, 1329:0, 31864:0, 10903:0,
       17325:h, 642:P*h, 9118:x7068_0-Q*P*h,
       9413:v[28730]//P, 28730:P*(v[28730]//P),
       21574:0, 1844:0}
t0=time.time(); ch,st=L.ripple(v,seeds)
print(f"ripple changed {len(ch)} vars ({time.time()-t0:.0f}s); x7068 preserved: {v[7068]==x7068_0}; "
      f"x7075={v[7075]} x21279={v[21279]} x2099==x9118: {v[2099]==v[9118]}")
B,F=rep(v,'after bit flip + local repair')
for a in B[:25]:
    print(f"   a{a} out={str(L.atom_out.get(a)):12s} in {len(L.atom2eq.get(a,{}))} eqs")
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','sw3_out.json'),'w'))
