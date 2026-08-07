"""2-bit messages around the best anchors, scored exactly in GF(p)."""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
anchors=[24601,2081,4287,13195]
t0=time.time(); out=[]
for A in anchors:
    res=[]
    for b in real:
        if b==A: continue
        v=msg({A,b}); F=fails(v)
        res.append((len(F),b,tree.get(b,'?')))
    res.sort()
    base=len(fails(msg({A})))
    print(f"anchor x{A} (alone: {base}):  best pairs -> "
          f"{[(b,n,tk) for n,b,tk in res[:6]]}   dist {dict(sorted(collections.Counter(n for n,_,_ in res).items())[:6])}"
          f"  ({time.time()-t0:.0f}s)", flush=True)
    out.append((A,res[:20]))
json.dump([[A,[[n,b,tk] for n,b,tk in r]] for A,r in out], open(os.path.join(HERE,'data','bits7.json'),'w'))
