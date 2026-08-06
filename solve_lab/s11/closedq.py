"""Rational diagnostic on the CLOSED region: is the p-factor still in the invariant
   once the move set is not artificially restricted?"""
import sys, os, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from ip9 import rational_solve
from closed import build_closed
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
for cap in [130, 300, 500, 900]:
    v=load_raw(src)
    v,FAIL,ROWS,used,M,rhs=build_closed(v,maxrows=cap,verbose=False)
    t0=time.time()
    x,piv=rational_solve(M,rhs)
    if x is None:
        print(f"  cap {cap:4d}: system {len(ROWS)}x{len(used)}  INCONSISTENT over Q ({time.time()-t0:.0f}s)", flush=True)
        continue
    D=1
    for t in x: D=D*t.denominator//math.gcd(D,t.denominator)
    pd = (D%P==0)
    cof = D//P if pd else D
    print(f"  cap {cap:4d}: system {len(ROWS)}x{len(used)}  D={len(str(D)):3d} digits  p|D={pd}  "
          f"cofactor={cof if cof<10**14 else str(cof)[:16]+'..'}  ({time.time()-t0:.0f}s)", flush=True)
