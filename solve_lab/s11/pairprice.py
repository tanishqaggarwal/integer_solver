import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
# the tight set (Hall violator) from the exhaustive matching, plus their atoms
TIGHT = {'a688':688,'a1618':1618,'a29539':29539,'a26731':26731,'a7881':7881,'a21050':21050,
         'a26839':26839,'a40065':40065,'a14445':14445,'a27139':27139,'a34580':34580,'a33796':33796}
MIRROR = [26719,26721,26723]
print("single-constraint equation costs:")
for k,a in sorted(TIGHT.items(), key=lambda kv: len(L.atom2eq.get(kv[1],{}))):
    print(f"  {k:8s} a{a}: {len(L.atom2eq.get(a,{}))} eqs")
mu=set()
for a in MIRROR: mu |= set(L.atom2eq.get(a,{}))
print(f"  mirror trio (26719,26721,26723): {len(mu)} eqs  <-- absorbs BOTH deficits")
print()
print("cheapest PAIRS by union of equations:")
res=[]
for (k1,a1),(k2,a2) in itertools.combinations(sorted(TIGHT.items()),2):
    u=set(L.atom2eq.get(a1,{})) | set(L.atom2eq.get(a2,{}))
    res.append((len(u),k1,k2))
res.sort()
for n,k1,k2 in res[:10]:
    print(f"  {k1}+{k2}: {n} eqs")
print(f"\n  mirror-trio option: {len(mu)} eqs")
