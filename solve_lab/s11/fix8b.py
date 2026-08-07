import sys, os, json, time, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
Q=7376877; LQ=15804267
CONST=[abs(c) for m,c in L.polys[3576].items() if m==(2081,)][0]
print("pin constant CONST == x6418 currently:", CONST==v[6418], " x26777 =",v[26777],
      " x26777 free:", 26777 not in L.definer, " x26777 atoms:",[a for a in L.var_atoms[26777]])
print("x2099 == x6418 :", v[2099]==v[6418])
x7068_0=v[7068]
R=(x7068_0-CONST)%LQ
g=math.gcd(Q*P%LQ, LQ)
print(f"need Q*P*h == x7068_0-CONST (mod {LQ});  gcd={g}; residue {R}; divisible: {R%g==0}")
if R%g==0:
    Mm=LQ//g
    h0=((R//g)*pow((Q*P//1)%LQ//g if False else (Q*P%LQ)//g, -1, Mm))%Mm
    # verify
    print("  check:", (Q*P*h0 - (x7068_0-CONST))%LQ==0)
    print("  h0 =",h0," M =",Mm)
