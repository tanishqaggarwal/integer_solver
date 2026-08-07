"""Are channel B's invariants degenerate, and can they vanish?

Within channel B (U=0, so no A/B bits) the four invariants depend on only 10-22 C bits and
10-16 D bits.  Sample that space and see how many distinct values they take and whether zero
appears -- the same question that made channel A's inv5 worth enumerating.
"""
import sys, os, json, time, random, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certsB.json')))]
dep=json.load(open(os.path.join(HERE,'data','chanB2.json')))
def inv(v): return tuple(sum(y*evalp(L.polys[a],v) for a,y in c.items())%P for c in CERT)
VB=[b for b in real if tree.get(b) in ('C','D')]
act=set()
for j in dep: act |= {b for b in dep[j] if tree.get(b) in ('C','D')}
act=sorted(act)
print(f"V-side bits moving any channel-B invariant: {len(act)} "
      f"(C {sum(1 for b in act if tree.get(b)=='C')}, D {sum(1 for b in act if tree.get(b)=='D')})")
rnd=random.Random(77)
vals=[collections.Counter() for _ in CERT]; zeros=collections.Counter()
t0=time.time(); N=0
best=None
for it in range(700):
    w=rnd.randrange(1,min(9,len(act))+1)
    S=set(rnd.sample(act,w))
    v=msg(S)
    if v[7715]!=0: continue          # must stay in channel B
    I=inv(v); N+=1
    nz=sum(1 for x in I if x==0); zeros[nz]+=1
    if best is None or nz>best[0]: best=(nz,sorted(S))
    for j,x in enumerate(I): vals[j][x]+=1
print(f"{N} channel-B messages ({time.time()-t0:.0f}s)")
for j in range(len(CERT)):
    print(f"  invB{j}: {len(vals[j])} distinct in {N}; zero seen {vals[j].get(0,0)}; "
          f"top multiplicity {vals[j].most_common(1)[0][1]}")
print("  #invariants simultaneously zero:", dict(sorted(zeros.items())))
print("  best message:", best)
