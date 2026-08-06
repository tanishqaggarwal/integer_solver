import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import atomval, load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
print("2458959 =", 2458959, "= 3 *", 2458959//3, " ; 819653 prime?",
      all(819653 % q for q in range(2, int(819653**0.5)+1)))
for nm,path in [('checkpoint 39026', os.path.join(HERE,'..','best','new_instance_partial_39026.json')),
                ('s11 best 39018',  os.path.join(HERE,'data','finish3_named.json'))]:
    v=load_raw(path)
    AV=[atomval(a,v) for a in range(L.NA)]
    F=[e for e in range(L.NEQ) if sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())!=0]
    vals=[sum(c*AV[a] for a,c in L.eq_atoms[e][2].items()) for e in F]
    print(f"\n{nm}: {len(F)} failing")
    print("   values divisible by P :", sum(1 for x in vals if x%P==0), "of", len(vals))
    print("   values divisible by 3 :", sum(1 for x in vals if x%3==0))
    g=0
    for x in vals: g=math.gcd(g,x)
    print("   gcd of failing values :", len(str(g)), "digits ; divisible by P:", g%P==0)
