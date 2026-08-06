import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
Q=7376877; LQ=15804267
CONST=[abs(c) for m,c in L.polys[3576].items() if m==(2081,)][0]
print("x3387 free:",3387 not in L.definer," atoms:",L.var_atoms[3387]," val:",v[3387])
print("x38744 == p:", v[38744]==P)
print()
print("x7068 %p  =", v[7068]%P)
print("CONST %p  =", CONST%P)
print("EQUAL:", v[7068]%P == CONST%P)
print("a22229 val %p =", L.evalpoly(L.polys[22229],v)%P)
print("(x7068-x2099) %p =", (v[7068]-v[2099])%P, " Q*x642 %p =", (Q*v[642])%P)
