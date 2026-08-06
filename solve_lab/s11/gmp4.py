"""Are the mirror checks FIXED POINTS or plain assignments?

a29539 says the free input x14853 must equal the computed x1308.  If x1308 does not depend on
x14853, this is a plain assignment.  If it does, it is a fixed-point equation -- and mod p an
affine fixed point is solved by one division, which is exactly the thing that has no integer
analogue.  Same question for a7930 (x24548 vs x25442) and for x9118/x8731.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
TGT=[7930,29539,35759,35760,40826,41512]
def state(u,val):
    v=list(base); v[u]=val%P; forwardp(v); return v
def resid(v): return [evalp(L.polys[a],v) for a in TGT]
def nfail(v): return sum(1 for a in CHK if evalp(L.polys[a],v))
print("base residues:", [str(x)[:14]+'..' for x in resid(base)], " failing:", nfail(base))
for u,partner in [(14853,1308),(24548,25442),(9118,None),(8731,None)]:
    r=[]
    for s in (0,1,2):
        v=state(u,(base[u]+s))
        r.append((v[partner] if partner else None, resid(v), nfail(v)))
    dep = (r[0][0]!=r[1][0]) if partner else None
    aff = all((r[2][1][j]-r[0][1][j])%P == 2*((r[1][1][j]-r[0][1][j])%P)%P for j in range(6))
    print(f"x{u}: partner x{partner} moves with it: {dep};  6-residue response affine: {aff};  "
          f"failing checks at +0/+1/+2 = {r[0][2]}/{r[1][2]}/{r[2][2]}")
    print("     d(residues)/du =", [str((r[1][1][j]-r[0][1][j])%P)[:12]+'..' for j in range(6)])
