"""Which bits move which invariant?

Invariant 5 came out identical for {x2081} and {x2081,x24601}, so the invariants are not all
functions of the whole message.  If each depends on only a few bits, hitting zero stops being a
2^256 problem and becomes several small ones.
"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
def inv(v): return tuple(sum(y*evalp(L.polys[a],v) for a,y in c.items())%P for c in CERT)
BASE={2081,24601}
v0=msg(BASE); I0=inv(v0)
print("base msg",sorted(BASE),"invariants nonzero:",[i for i,x in enumerate(I0) if x])
dep=collections.defaultdict(set)
t0=time.time()
for b in real:
    S=(BASE-{b}) if b in BASE else (BASE|{b})
    if not S: continue
    v=msg(S); I=inv(v)
    for j in range(len(CERT)):
        if I[j]!=I0[j]: dep[j].add(b)
print(f"scanned ({time.time()-t0:.0f}s)")
for j in range(len(CERT)):
    d=sorted(dep[j])
    byt=collections.Counter(tree.get(b,'?') for b in d)
    print(f"  invariant {j} (support {len(CERT[j])}): moved by {len(d)} bits  by tree {dict(byt)}")
    if len(d)<=45: print(f"      bits: {d}")
json.dump({str(j):sorted(dep[j]) for j in dep}, open(os.path.join(HERE,'data','bits11.json'),'w'))
