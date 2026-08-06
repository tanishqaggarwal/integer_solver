"""Clear the obstruction WITHOUT disturbing the copy network.

fix2 moved x642 to a multiple of p, which moved x7068 = x2099 + 7376877*x642, and x7068 is
copied into free 'mirror' inputs all over the circuit -- hence the cascade.

But x2099 = x6418 (the other two summands vanish: x10878 = x6788*x31861 = 0 and
x25297 = x9118*x21279 = 0), and x6418 is a FREE load-pin input whose pin

    a3576 :  x2081*(x6418 - HUGE) = 15804267*x26777        (x2081 = 1)

quantises it modulo 15804267 -- NOT modulo p.  So choose

    x642 = p*h   with   7376877*(p*h - x642_orig) == 0  (mod 15804267)
    x6418 = x7068_orig - 7376877*p*h

and x7068 never moves.  The copy network is untouched.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v=load_raw(src)
def rep(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken atoms={len(B)} {B[:20]}  failing={len(F)} score={L.NEQ-len(F)}")
    for a in B[:20]: print(f"      a{a} out={str(L.atom_out.get(a)):12s} in {len(L.atom2eq.get(a,{}))} eqs")
    return B,F
rep(v,'start')
Q=7376877; LQ=15804267
x642_0=v[642]; x7068_0=v[7068]; HUGE=v[6418]
d=__import__('math').gcd(Q,LQ); M=LQ//d
print(f"gcd(7376877,15804267)={d}  modulus M={M}")
h0=(x642_0 * pow(P,-1,M)) % M
t=(x642_0//P - h0)//M
h=h0+M*t
assert (Q*(P*h-x642_0))%LQ==0, "pin congruence"
print(f"h chosen; p*h has {len(str(abs(P*h)))} digits vs x642 {len(str(abs(x642_0)))}")
new6418 = x7068_0 - Q*P*h
assert (new6418-HUGE)%LQ==0, "load pin divisibility"
seeds={
  17325:h, 642:P*h,
  6418:new6418, 26777:(new6418-HUGE)//LQ,
  9118: v[9118]-(v[9118]%P), 8731: v[8731]-(v[8731]%P),
  21574:0, 1844:0,
}
seeds[1329]  = 5113045*v[7075]*seeds[9118]//P
seeds[29854] = 5113045*v[7075]*seeds[9118]
seeds[31864] = -v[7075]*seeds[8731]
seeds[10903] = seeds[31864]//P
seeds[9413]  = v[28730]//P
seeds[28730] = P*seeds[9413]
t0=time.time(); ch,st=L.ripple(v,seeds)
print(f"ripple changed {len(ch)} vars ({time.time()-t0:.0f}s);  x7068 unchanged: {v[7068]==x7068_0}")
B,F=rep(v,'after')
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','fix8_out.json'),'w'))
