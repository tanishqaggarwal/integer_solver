"""What kind of atoms make up the obstruction certificate?

A load pin reads  bit*(x - C) - m*handle.  With the bit OFF the pin degenerates to
handle = 0 -- it LOCKS the handle.  With the bit ON it reads x = C + m*handle, which frees the
handle to be anything.  So a certificate containing pins of currently-OFF bits is telling us
which bits to switch on to unlock knobs.
"""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
SUP=[42245,36040,29253,17518,8948,3580,8824,14515,27069,14523,1645,29457,21617,15462,31643,
     36289,18593,15753,17931,33539,28434,7289,3254,7281,1637,15456,13198,29824,17923,1653,
     7930,3568,3578,3584,14507,3260,25676,33792,40623,40562,35856,1739,29165,1524,17524]
pins=collections.defaultdict(list)
for a in SUP:
    Pp=L.polys[a]
    bs=[u for u in L.avars[a] if u in BITS]
    # pin shape: some monomial (bit, x) and a constant*bit and a small*handle
    shape=None
    for m,c in Pp.items():
        if len(m)==1 and m[0] in BITS and abs(c)>10**60: shape=(m[0],c)
    kind = 'PIN' if shape else ('bit-check' if len(Pp)==2 and any(len(m)==2 for m in Pp) and bs else 'other')
    fv=[u for u in sorted(L.avars[a]) if u in FREE]
    print(f"a{a:6d} {kind:10s} nterms={len(Pp):3d} bits={bs} bitvals={[base[u]%P for u in bs]} "
          f"free={fv[:6]}{'..' if len(fv)>6 else ''}")
    if shape: pins[shape[0]].append(a)
print()
print("bits owning pins in the certificate:", {u:(base[u]%P, v) for u,v in pins.items()})
