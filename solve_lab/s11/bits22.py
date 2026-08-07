"""Systematic enumeration of invariant 5 over the 18 bits that move it."""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
ACT=[1530,2081,4287,6821,11368,12054,13195,14808,16586,17406,17760,21074,22351,22562,23751,24365,34974,35979]
INERT=[438,490,1603,2455,5643,5910,7401,12095,14823,16348,16827,18022,22623,24517,27393,28005,33287,37862,38625]
t0=time.time(); vals=collections.Counter(); zero=[]
n=0
for w in range(0,6):
    for S in itertools.combinations(ACT,w):
        x=inv5(msg({24601}|set(S))); vals[x]+=1; n+=1
        if x==0: zero.append(S); print("  *** inv5 == 0 at", S, flush=True)
    print(f"  weight {w} done: {n} subsets, {len(vals)} distinct, zero={len(zero)} ({time.time()-t0:.0f}s)", flush=True)
print(f"TOTAL {n} subsets of the 18 active bits (weight <= 5): {len(vals)} distinct inv5 values, zero found {len(zero)}")
print("  multiplicity profile:", collections.Counter(vals.values()).most_common(8))
json.dump({'n':n,'distinct':len(vals),'zero':[list(z) for z in zero]},
          open(os.path.join(HERE,'data','bits22.json'),'w'))
