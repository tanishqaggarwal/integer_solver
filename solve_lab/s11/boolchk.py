import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
LD=json.load(open(os.path.join(HERE,'data','loads.json')))['loads']
BITS=sorted(int(b) for b in LD)
# look for atoms enforcing b^2 = b  (or b^2 - b, b*(b-1))
found=collections.defaultdict(list)
for a in range(L.NA):
    Pp=L.polys[a]
    for mm,c in Pp.items():
        if len(mm)==2 and mm[0]==mm[1] and mm[0] in set(BITS):
            found[mm[0]].append(a)
print(f"atoms containing bit^2 for a message bit: {len(found)} bits affected")
for b in list(found)[:6]:
    print(f"  x{b}: atoms {found[b][:5]}")
# check a sample bit directly: does setting it to 2 break anything beyond its own pins?
import copy
for b in BITS[:4]:
    v=[0]*L.NVARS
    fw.forward(v)
    base=set(fw.bad_checks(v))
    v[b]=2; fw.forward(v)
    nb=[a for a in fw.bad_checks(v) if a not in base]
    v[b]=1; fw.forward(v)
    nb1=[a for a in fw.bad_checks(v) if a not in base]
    print(f"  bit x{b}: value 2 -> {len(nb)} new bad ; value 1 -> {len(nb1)} new bad")
