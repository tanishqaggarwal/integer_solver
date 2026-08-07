"""Invariant 5 sees only 18 C bits.  Is it additive on them?

If yes, all 2^18 = 262,144 subsets can be evaluated arithmetically -- no forward evaluations at
all -- and we can say definitively whether inv5 = 0 is reachable.
"""
import sys, os, json, time, itertools, random, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
ACT=[1530,2081,4287,6821,11368,12054,13195,14808,16586,17406,17760,21074,22351,22562,23751,24365,34974,35979]
print("active C bits for inv5:",len(ACT), flush=True)
i0=inv5(msg({24601}))
d={b:(inv5(msg({24601,b}))-i0)%P for b in ACT}
rnd=random.Random(31); ok=bad=0; mism=[]
for _ in range(40):
    k=rnd.randrange(2,6)
    S=rnd.sample(ACT,k)
    got=inv5(msg({24601}|set(S)))
    pred=(i0+sum(d[b] for b in S))%P
    if got==pred: ok+=1
    else: bad+=1; mism.append((tuple(S),got,pred))
print(f"additive on {ok}/{ok+bad} random subsets of the 18 active bits", flush=True)
if bad==0:
    print("=> ADDITIVE.  Enumerating all 2^18 subset sums arithmetically...", flush=True)
    t0=time.time()
    vals=[i0]
    for b in ACT:
        db=d[b]
        vals=vals+[(x+db)%P for x in vals]
    st=set(vals)
    print(f"   {len(vals)} sums in {time.time()-t0:.1f}s;  distinct {len(st)};  ZERO present: {0 in st}")
    if 0 in st:
        for m in range(1<<len(ACT)):
            s=i0
            for i,b in enumerate(ACT):
                if m>>i&1: s=(s+d[b])%P
            if s==0:
                print("   *** inv5 = 0 for C-subset", [ACT[i] for i in range(len(ACT)) if m>>i&1]); break
else:
    print("  example mismatches:", [(S,str(g)[:12],str(p)[:12]) for S,g,p in mism[:3]])
