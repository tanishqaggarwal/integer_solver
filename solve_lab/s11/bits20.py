"""What statistic of the C-subset does invariant 5 see?

37 single-bit C-subsets give only ~15 distinct values, so the bits fall into equivalence classes.
If inv5 sees only a coarse statistic of the subset, the search collapses further.
"""
import sys, os, json, time, itertools, collections, random
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
Cb=sorted(b for b in real if tree.get(b)=='C')
single={b:inv5(msg({24601,b})) for b in Cb}
cls=collections.defaultdict(list)
for b,x in single.items(): cls[x].append(b)
print(f"{len(cls)} classes among the 37 C bits under inv5 (weight-1):")
for x,bs in sorted(cls.items(), key=lambda z:-len(z[1])):
    print(f"   {len(bs):2d} bits -> {str(x)[:18]}..   {bs}")
# does inv5(S) depend only on the multiset of classes?
rnd=random.Random(17)
byclass={x:bs for x,bs in cls.items()}
tests=0; agree=0
for _ in range(40):
    k=rnd.randrange(2,5)
    xs=rnd.sample(list(byclass),k)
    S1={rnd.choice(byclass[x]) for x in xs}
    S2={rnd.choice(byclass[x]) for x in xs}
    if S1==S2: continue
    a=inv5(msg({24601}|S1)); b=inv5(msg({24601}|S2))
    tests+=1; agree+= (a==b)
print(f"\nsame class-multiset => same inv5:  {agree}/{tests}")
