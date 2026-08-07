"""Is invariant 5 additive over the 37 C-tree bits?

It depends on those 37 and nothing else.  If it is additive in them, finding a zero is a
subset-sum with 37 items in GF(p) -- meet-in-the-middle over 2^18/2^19 lists, entirely feasible.
"""
import sys, os, json, time, itertools, random, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
Cbits=sorted(json.load(open(os.path.join(HERE,'data','bits11.json')))['5'])
print(f"C-tree bits: {len(Cbits)}")
ANCH={24601}
def I(S): return inv5(msg(ANCH|set(S)))
i_none=I([])                      # V=0 here; still well defined
singles={b:I([b]) for b in Cbits}
print("inv5 with no C bit (V=0):", str(i_none)[:20]+'..')
print("inv5 for single C bits (first 5):", [str(singles[b])[:14]+'..' for b in Cbits[:5]])
ok=bad=0; t0=time.time()
rnd=random.Random(7)
for _ in range(30):
    k=rnd.randrange(2,5)
    S=rnd.sample(Cbits,k)
    got=I(S)
    pred=(sum(singles[b] for b in S)-(len(S)-1)*i_none)%P
    if got==pred: ok+=1
    else: bad+=1
print(f"additive on {ok}/{ok+bad} random subsets of size 2-4 ({time.time()-t0:.0f}s)")
if bad:
    S=rnd.sample(Cbits,2)
    print("  example mismatch check, pair",S, str(I(S))[:16]+'..', "vs pred",
          str((singles[S[0]]+singles[S[1]]-i_none)%P)[:16]+'..')
