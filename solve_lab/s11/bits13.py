"""How many distinct values do the six invariants take over random messages?

If they behave like random elements of GF(p), hitting zero is hopeless by search.  If they take
few values, the search collapses.  This is the question that decides whether the discrete route
is worth pursuing at all.
"""
import sys, os, json, time, random, collections
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
rnd=random.Random(11)
vals=[collections.Counter() for _ in CERT]
zeros=collections.Counter()
t0=time.time(); N=0
best=None
for it in range(500):
    w=rnd.choice([1,2,2,3,3,4,6,10,20,40])
    S=set(rnd.sample(real,w))
    # keep the OR gate satisfied
    if not (set(S)&set(byt['A']+byt['B'])) and not (set(S)&set(byt['C']+byt['D'])): continue
    v=msg(S); I=inv(v); N+=1
    nz=sum(1 for x in I if x==0)
    zeros[nz]+=1
    if best is None or nz>best[0]: best=(nz,sorted(S))
    for j,x in enumerate(I): vals[j][x]+=1
print(f"{N} random messages ({time.time()-t0:.0f}s)")
for j in range(len(CERT)):
    print(f"  invariant {j}: {len(vals[j])} distinct values in {N} samples; "
          f"most common count {vals[j].most_common(1)[0][1]}; zero seen: {vals[j].get(0,0)}")
print("  #invariants simultaneously zero:", dict(sorted(zeros.items())))
print("  best message so far:", best)
