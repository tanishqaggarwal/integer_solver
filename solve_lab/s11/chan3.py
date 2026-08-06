"""Reach the third MUX channel, x34606 = U*(1-V).

x7715 = OR(x8599,x21839) = U ; x34554 = OR(x7304,x25956) = V
x15298 = U*V   x5647 = (1-U)*V   x34606 = U*(1-V)
The checkpoint sits at U=V=1 and the 39,018 state at U=0,V=1.  U=1,V=0 has never been built.
Turning off whichever ON bit feeds V should land there.
"""
import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from ip7 import load_raw
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000); LAB=os.path.join(HERE,'..')
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
def fails(v): return [a for a in CHK if evalp(L.polys[a],v)]
def eqfail(v):
    AV=[evalp(L.polys[a],v) for a in range(L.NA)]
    return [e for e in range(L.NEQ) if sum(c*AV[a] for a,c in L.eq_atoms[e][2].items())%P]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
def rep(v,tag):
    F=fails(v); E=eqfail(v)
    print(f"{tag}: U={v[7715]} V={v[34554]} | x15298={v[15298]} x5647={v[5647]} x34606={v[34606]} "
          f"| checks fail {len(F)} {F[:10]} | eqs fail {len(E)} ceiling {L.NEQ-len(E)}")
    return F
rep(base,'base (U=1,V=1)')
for b in (2081,24601):
    v=list(base); v[b]=0; forwardp(v)
    rep(v,f'bit x{b} OFF')
v=list(base); v[2081]=0; v[24601]=0; forwardp(v)
rep(v,'both OFF')
