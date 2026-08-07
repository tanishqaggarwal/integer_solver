"""Search the C-tree subsets for invariant 5 = 0.

inv5 is moved only by C-tree bits, at every anchor tested -- so it is a function of the C-subset
alone (37 bits), independent of A, B and D.  That splits the search: find a C-subset zeroing
inv5, then use the other 219 bits on the remaining invariants.
"""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
Cb=sorted(b for b in real if tree.get(b)=='C')
ANCH={24601}
print(f"C bits: {len(Cb)}")
t0=time.time(); vals=collections.Counter(); zero=[]
allS=[()]+[(b,) for b in Cb]+list(itertools.combinations(Cb,2))+list(itertools.combinations(Cb,3))
print(f"enumerating {len(allS)} C-subsets of weight <= 3")
for i,S in enumerate(allS):
    v=msg(ANCH|set(S)); x=inv5(v)
    vals[x]+=1
    if x==0: zero.append(S); print(f"  *** inv5 == 0 at C-subset {S}", flush=True)
    if i%1500==0: print(f"   {i}/{len(allS)} distinct={len(vals)} ({time.time()-t0:.0f}s)", flush=True)
print(f"done ({time.time()-t0:.0f}s): {len(vals)} distinct inv5 values over {len(allS)} subsets")
print("  most common values and their multiplicity:", [(str(k)[:12]+'..',c) for k,c in vals.most_common(6)])
print("  zero found:", len(zero), zero[:5])
json.dump({'zero':[list(z) for z in zero],'ndistinct':len(vals)}, open(os.path.join(HERE,'data','bits16.json'),'w'))
