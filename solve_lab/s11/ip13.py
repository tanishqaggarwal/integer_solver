"""IP #13 -- confirm the TRUE minimal invariant at the state where D/p = 1."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip8 import build
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
src=os.path.join(HERE,'data','closehit2.json')
v=load_raw(src)
print("=== closehit2.json (D/p == 1 state)")
v,FAIL,used,M,rhs,nf=build(v)
for d,lab in [(1,'1'),(2,'2'),(3,'3'),(P,'P')]:
    t0=time.time()
    x=solve_int(M,[d*r for r in rhs])
    print(f"   M x = {lab:2s} * rhs : {'SOLVABLE' if x else 'no':9s} ({time.time()-t0:.0f}s)", flush=True)
    if x:
        print(f"   => the ENTIRE obstruction at this state is a single divisibility by {lab}")
        break
