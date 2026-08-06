"""Turn x2081 OFF to unlock x6418 and x12553.

Non-bit movers of the four failing checks:
    a7930  <- x12553, x24548          a29539 <- x6418, x14853
    a40826 <- x6418, x14515, x14853, x19750    a41512 <- x96, x9280, x12553, x18027, x24548, x27711
x6418 and x12553 are the two values LOADED by bit x2081.  With the bit ON their pins
(a3576, a3578) force them to fixed residues mod p, so they are frozen.  With the bit OFF the
pins degenerate to handle == 0 and both become free mod p -- and between them they move all
four failing checks.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
def fails(v): return [a for a in CHK if evalp(L.polys[a],v)]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
print("base (x2081=1):", len(fails(base)), fails(base))
v=list(base); v[2081]=0; forwardp(v)
F=fails(v)
print(f"x2081 -> 0: failing = {len(F)}  {F[:20]}")
print("  OR-gate atom a23000 =", evalp(L.polys[23000],v), " (must be 0)")
# are x6418 / x12553 free mod p now?
for u in (6418,12553):
    w=list(v); w[u]=(w[u]+1)%P; forwardp(w)
    d=[a for a in CHK if evalp(L.polys[a],w)!=evalp(L.polys[a],v)]
    print(f"  x{u} now moves checks: {d}")
json.dump([int(x) for x in v], open(os.path.join(HERE,'data','gmp22_off.json'),'w'))
