import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
def rep(v,tag):
    AV=[L.evalpoly(L.polys[a],v) for a in range(L.NA)]
    B=[a for a in range(L.NA) if AV[a]!=0]; F=L.failing_eqs(AV)
    print(f"{tag}: broken={B} failing={len(F)} score={L.NEQ-len(F)}")
    for a in B: print(f"    a{a} in {len(L.atom2eq.get(a,{}))} eqs")
    return B,F
v=load_raw(os.path.join(HERE,'data','modp3_out.json'))
print("x15298=",v[15298]," x5101==p:",v[5101]==P," x11150%p==0:",v[11150]%P==0)
seeds={}
n=-v[11150]*v[15298]                      # want x4007 = n = p*x30317
print("  a19297: p | needed x4007 :", n%P==0)
if n%P==0: seeds[30317]=n//P
n=v[15298]*v[25739]                        # want 6672769*x29804 = n, x29804 = p*x5146
print("  a19299: 6672769*p | needed :", n%(6672769*P)==0)
if n%(6672769*P)==0: seeds[5146]=n//(6672769*P)
n=537773*v[15298]*v[37758]                 # want x35605 = n = p*x2936
print("  a30984: p | needed x35605 :", n%P==0)
if n%P==0: seeds[2936]=n//P
L.ripple(v,seeds)
B,F=rep(v,'after correct handles')
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','modp5_out.json'),'w'))
