import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine, deep
from zsolve import solve_int
BAD=[26719,26721,26723,26733,28438,32342,36185]
v=engine.apply_theta({5096:0,19750:0})
locked=set(engine.BITS)|engine.DERIVED
cands=set([5096,19750])
for a in BAD:
    h,base=deep.handles(v,a,locked=locked)
    print(f"a{a}: handles={[(t,len(L.var_atoms[t])) for t,_ in h]}")
    for t,_ in h: cands.add(t)
cands-=engine.DERIVED
print("controls:", sorted(cands))
json.dump(sorted(cands), open('controls.json','w'))
