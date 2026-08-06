"""Exhaustive weight-2 message scan: all C(256,2) = 32,640 patterns, scored exactly in GF(p)."""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
t0=time.time(); res=[]
for i,(a,b) in enumerate(itertools.combinations(real,2)):
    v=msg({a,b}); F=fails(v)
    res.append((len(F),a,b,tuple(F) if len(F)<=8 else ()))
    if i%4000==0:
        print(f"   {i}/32640 best={min(r[0] for r in res)} ({time.time()-t0:.0f}s)", flush=True)
res.sort()
print(f"done ({time.time()-t0:.0f}s).  distribution: {dict(sorted(collections.Counter(r[0] for r in res).items())[:12])}")
print("best weight-2 messages:")
for n,a,b,F in res[:25]:
    print(f"   {{x{a}, x{b}}} trees {tree.get(a,'?')}{tree.get(b,'?')}  failing = {n}  {list(F)}")
json.dump([[r[0],r[1],r[2]] for r in res[:400]], open(os.path.join(HERE,'data','bits8.json'),'w'))
