import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
GATES=[a for a in range(L.NA) if L.atom_out.get(a) is not None]
cost={a:len(L.atom2eq.get(a,{})) for a in GATES}
cheap=sorted(GATES,key=lambda a:cost[a])[:60]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
bd={a:evalp(L.polys[a],base) for a in CHK}
print("cheapest gate atoms:", [(a,cost[a]) for a in cheap[:14]])
for g in cheap[:14]:
    t=L.atom_out[g][1]
    v=list(base); v[t]=(v[t]+1)%P; forwardp_frozen(v,{t})
    d=[a for a in CHK if evalp(L.polys[a],v)!=bd[a]]
    print(f"  a{g} cost={cost[g]} frees x{t}: moves {len(d)} checks {d[:12]}")
