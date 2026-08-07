import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
for u in [27676,7574,31861,14865,9413,28730,4432,19964]:
    print(f"x{u} {'FREE' if u not in L.definer else 'DEF a%d'%L.definer[u]} val={sh(v[u])} atoms={len(L.var_atoms[u])}")
print()
for a in [3568,3570]:
    for m,c in L.polys[a].items():
        print(f"  a{a}: {c} * {'*'.join('x%d'%u for u in m) if m else '1'}")
    print()
# what does x7075=0 / x21279=1 imply?  simulate the pure bit flip with pins satisfied
import copy
BIG1=[abs(c) for m,c in L.polys[3568].items() if m==(4287,)][0]
BIG2=[abs(c) for m,c in L.polys[3570].items() if m==(4287,)][0]
print("BIG1 (pin for x31861) digits:",len(str(BIG1)), " BIG1 % 13479571 =", BIG1%13479571)
print("BIG2 (pin for x14865) digits:",len(str(BIG2)))
