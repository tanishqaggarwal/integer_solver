"""Characterise the invariant value sets.

500 random messages gave only 88-176 distinct values per invariant instead of 500 -- so the
invariants live on a small discrete set, not on all of GF(p).  Map that set: does it saturate,
are the invariants proportional to each other, and is 0 in it?
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
UB=byt['A']+byt['B']; VB=byt['C']+byt['D']
rnd=random.Random(23)
vals=[dict() for _ in CERT]      # value -> example message
t0=time.time(); N=0
seen=set()
for it in range(4000):
    w=rnd.choice([1,2,2,3,3,4,5,6,8,12,20,40,80])
    S=frozenset(rnd.sample(real,w))
    if S in seen: continue
    if not (S&set(UB)) and not (S&set(VB)): continue
    seen.add(S)
    v=msg(S); I=inv(v); N+=1
    for j,x in enumerate(I):
        if x not in vals[j]: vals[j][x]=sorted(S)
    if N%800==0: print(f"   {N} msgs, distinct per invariant {[len(d) for d in vals]} ({time.time()-t0:.0f}s)", flush=True)
print(f"{N} messages ({time.time()-t0:.0f}s)")
for j in range(len(CERT)):
    print(f"  invariant {j}: {len(vals[j])} distinct values; zero present: {0 in vals[j]}")
# proportionality between invariants
import itertools
pairs=[]
sample=[]
rnd2=random.Random(5)
for S in list(seen)[:60]:
    v=msg(set(S)); sample.append(inv(v))
for i,j in itertools.combinations(range(len(CERT)),2):
    rs={ (a[i]*pow(a[j],-1,P))%P for a in sample if a[j] }
    if len(rs)==1: pairs.append((i,j,rs.pop()))
print("  proportional invariant pairs (i,j,ratio):", [(i,j,str(r)[:14]+'..') for i,j,r in pairs])
json.dump({str(j):[str(x) for x in sorted(vals[j])[:500]] for j in range(len(CERT))},
          open(os.path.join(HERE,'data','bits14_vals.json'),'w'))
