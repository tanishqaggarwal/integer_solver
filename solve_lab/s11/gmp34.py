"""Additivity WITHIN a fixed channel.

The pairwise test failed only on x1308, whose dependence runs through x3896 = x7304*x25956 and
the OR-trees -- i.e. through the channel/mirror INDICATORS, which are boolean functions of which
trees have a bit on.  Restricted to bits that leave U, V, x3896 and x38170 unchanged, the map
should be affine in the loaded constants -- and then the bit search is a subset-sum, not a
2^256 walk.
"""
import sys, os, json, time, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
PIN=collections.defaultdict(list)
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=3: continue
    bb=[m[0] for m,c in Pp.items() if len(m)==1 and m[0] in BITS and abs(c)>10**60]
    if not bb: continue
    b=bb[0]
    q=[m for m in Pp if len(m)==2 and b in m]
    if not q: continue
    x=[t for t in q[0] if t!=b][0] if q[0][0]!=q[0][1] else q[0][0]
    C=-[c for m,c in Pp.items() if m==(b,)][0]
    hs=[m[0] for m,c in Pp.items() if len(m)==1 and m[0]!=b]
    PIN[b].append((x,C%P,hs[0] if hs else None))
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
IND=[7715,34554,3896,38170,15298,5647,34606]
i0=[base[u]%P for u in IND]
def turn_on(v,b):
    v[b]=1
    for x,C,h in PIN.get(b,()):
        v[x]=C
        if h in FREE: v[h]=0
    forwardp(v); return v
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
same=[]
for b in real:
    if base[b]%P==1: continue
    v=list(base); turn_on(v,b)
    if [v[u]%P for u in IND]==i0: same.append(b)
print(f"bits that leave the channel/mirror indicators unchanged: {len(same)} of {len(real)}")
PROBE=[25442,1308,27522,11150,25739,37758,6858,25295,12000,12926,21364]
v0=[base[u]%P for u in PROBE]
single={}
for b in same[:14]:
    v=list(base); turn_on(v,b); single[b]=[v[u]%P for u in PROBE]
ok=bad=0
for i in range(min(8,len(single))):
    for j in range(i+1,min(9,len(single))):
        bi,bj=list(single)[i],list(single)[j]
        v=list(base); turn_on(v,bi); turn_on(v,bj)
        vv=[v[u]%P for u in PROBE]
        pred=[(single[bi][k]+single[bj][k]-v0[k])%P for k in range(len(PROBE))]
        mm=[PROBE[k] for k in range(len(PROBE)) if vv[k]!=pred[k]]
        if mm: bad+=1
        else: ok+=1
        print(f"   x{bi}+x{bj}: {'ADDITIVE' if not mm else 'mismatch at '+str(mm)}")
print(f"additive pairs: {ok}, non-additive: {bad}")
