"""Who can move each critical check, mod p?"""
import sys, os, json, pickle, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
cols=D['cols']; base=D['base']; bd=D['bd']
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
for a in [7930,29539,40826,41512,19297,19299,30984,21617,33796,36185,40812,37662,25676,42245]:
    mv=[u for u,d in cols.items() if a in d]
    nb=[u for u in mv if u not in BITS]
    print(f"a{a:6d} {'FAIL' if bd.get(a) else ' ok '}: {len(mv):3d} movers, {len(nb):3d} non-bit -> {nb[:12]}")
