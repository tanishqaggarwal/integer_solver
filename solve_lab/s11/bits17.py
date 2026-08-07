"""Image of invariant 5 over C-subsets of ALL weights -- does it grow, and does it contain 0?"""
import sys, os, json, time, random, collections
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
rnd=random.Random(99)
vals=collections.Counter(); t0=time.time(); zero=[]
byw=collections.defaultdict(set)
for it in range(1600):
    w=rnd.randrange(0,len(Cb)+1)
    S=set(rnd.sample(Cb,w))
    v=msg({24601}|S); x=inv5(v)
    vals[x]+=1; byw[w].add(x)
    if x==0: zero.append(sorted(S)); print("  *** inv5 == 0 at",sorted(S), flush=True)
    if it%400==0: print(f"   {it} distinct={len(vals)} ({time.time()-t0:.0f}s)", flush=True)
print(f"1600 random C-subsets of all weights: {len(vals)} distinct inv5 values ({time.time()-t0:.0f}s)")
print("  top multiplicities:", [c for _,c in vals.most_common(5)])
print("  distinct values by subset weight (sample):",
      {w:len(byw[w]) for w in sorted(byw)[:12]})
print("  zero found:", len(zero))
