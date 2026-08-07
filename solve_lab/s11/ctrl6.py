import sys, os, json, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, sys6
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
LD=json.load(open(os.path.join(HERE,'data','loads.json')))['loads']
BITSET=set(int(b) for b in LD)
S=json.load(open(os.path.join(HERE,'data','sys6.json')))
NAMES=sys6.NAMES
POOL=sorted({u for nm in NAMES for u in S[nm] if u not in BITSET})
rnd=random.Random(3)
th={u: sys6.BASE[u]+rnd.randrange(1,1<<40) for u in POOL}
r=sys6.six(sys6.ev(th))
for i,nm in enumerate(NAMES):
    nz=[]
    for c in POOL:
        t2=dict(th); t2[c]=th[c]+1
        if (sys6.six(sys6.ev(t2))[i]-r[i])%P: nz.append(c)
    print(f"{nm:8s}: {nz}")
