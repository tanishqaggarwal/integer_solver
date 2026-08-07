"""True dependency: from a clean anchor, which bits actually move each invariant?"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
def inv(v): return tuple(sum(y*evalp(L.polys[a],v) for a,y in c.items())%P for c in CERT)
byt=collections.defaultdict(list)
for b in real: byt[tree.get(b,'?')].append(b)
for ANCH,name in [({24601},'{x24601}  (U=1,V=0, channel C)'),
                  ({2081},'{x2081}   (U=0,V=1, channel B)'),
                  ({2081,24601},'{x2081,x24601} (channel A)')]:
    I0=inv(msg(ANCH))
    dep=collections.defaultdict(list)
    for b in real:
        S=(ANCH-{b}) if b in ANCH else (ANCH|{b})
        if not S: continue
        I=inv(msg(S))
        for j in range(len(CERT)):
            if I[j]!=I0[j]: dep[j].append(b)
    print(f"anchor {name}")
    for j in range(len(CERT)):
        d=dep[j]; byq=collections.Counter(tree.get(b,'?') for b in d)
        print(f"   inv{j}: moved by {len(d):3d} bits {dict(byq)}")
    print()
