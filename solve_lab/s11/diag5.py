import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
KEY=[35756,35757,35758,35759,35760,35761,35762,22229,22230]
print("equations touching these atoms:")
allq=set()
for a in KEY: allq |= set(L.atom2eq.get(a,{}))
print(sorted(allq))
AV=[atomval(a,v) for a in range(L.NA)]
def eqv(e): return sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())
for e in sorted(allq):
    m,sq,co=L.eq_atoms[e]
    print(f"  eq{e} sq={sq} val={'FAIL' if eqv(e)!=0 else 'ok'} coeffs on KEY: "
          f"{ {a:co[a] for a in KEY if a in co} }  others:{ {a:c for a,c in co.items() if a not in KEY} }")
print()
A = v[7068]-v[2099]
B = v[28730]
S = v[1956]*v[17065]
T = 5113045*v[7075]*v[9118]
U = v[7075]*v[8731]
for nm,x in [('A',A),('B',B),('S',S),('T',T),('U',U)]:
    print(f"{nm}: digits={len(str(abs(x)))}  mod p = {x%P}   ==0 mod p? {x%P==0}")
print("A mod 7376877 =", A%7376877, " A/7376877 exact?", A%7376877==0)
print("x28730 % P =", v[28730]%P)
