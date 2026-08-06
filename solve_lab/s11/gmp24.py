"""Are the message bits ADDITIVE on the residues?

Each ON bit b pins its loaded variable to a fixed constant C_b (mod p) and switches its path on.
If the quantities the mirror checks compare depend additively on which bits are on, then finding
the right pattern is a subset-sum in GF(p) rather than a blind 2^256 search.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
CHK=[a for a in range(L.NA) if L.atom_out.get(a) is None]
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
# pins:  b*x - C*b - m*h
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
    hs=[(m[0],c) for m,c in Pp.items() if len(m)==1 and m[0]!=b]
    PIN[b].append((a,x,C,hs[0] if hs else None))
print("bits with a pin:",len(PIN))
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
forwardp(base)
def turn_on(v,b):
    v[b]=1
    for a,x,C,h in PIN.get(b,()):
        v[x]=C%P
        if h and h[0] in FREE: v[h[0]]=0
    forwardp(v)
    return v
PROBE=[25442,1308,27522,11150,25739,37758,6858,25295]
def vals(v): return [v[u]%P for u in PROBE]
v0=vals(base)
cands=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))][:6]
print("probe vars:",PROBE)
print("base:",[str(x)[:10]+'..' for x in v0])
single={}
for b in cands:
    v=list(base); turn_on(v,b); single[b]=vals(v)
    print(f"  +x{b}: delta nonzero at {[i for i in range(len(PROBE)) if (single[b][i]-v0[i])%P]}")
print("\nadditivity test on pairs:")
for i in range(len(cands)):
    for j in range(i+1,len(cands)):
        bi,bj=cands[i],cands[j]
        v=list(base); turn_on(v,bi); turn_on(v,bj)
        vv=vals(v)
        pred=[(single[bi][k]+single[bj][k]-v0[k])%P for k in range(len(PROBE))]
        ok=all(vv[k]==pred[k] for k in range(len(PROBE)))
        print(f"   x{bi}+x{bj}: additive = {ok}   mismatched components {[k for k in range(len(PROBE)) if vv[k]!=pred[k]]}")
