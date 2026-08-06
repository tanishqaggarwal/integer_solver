"""Correct test: verify the certificates at the EXACT state where they were derived,
then at a sibling message reached by the minimal change (swap the C bit and its two loads)."""
import sys, os, json, time, random, pickle
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp26 import forwardp_frozen
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
D=json.load(open(os.path.join(HERE,'data','bits1.json')))
PIN={int(k):[(a,x,int(C)) for a,x,C in v] for k,v in D['pins'].items()}
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
R=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
BITS=set()
for a in range(L.NA):
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    ks=sorted(Pp.keys(), key=len)
    if len(ks[0])==1 and len(ks[1])==2 and ks[1]==(ks[0][0],ks[0][0]) and Pp[ks[0]]==-Pp[ks[1]]:
        BITS.add(ks[0][0])
live=[u for u in R['cols'] if u not in BITS]
rows=set()
for c in CERT: rows |= set(c)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
def swap(v, old, new):
    v=list(v); v[old]=0; v[new]=1
    for a,x,C in PIN[old]: v[x]=0
    for a,x,C in PIN[new]:
        v[x]=C%P
        Pp=L.polys[a]
        hs=[m[0] for m,c in Pp.items() if len(m)==1 and m[0]!=new]
        if hs and hs[0] in FREE: v[hs[0]]=0
    forwardp(v); return v
rnd=random.Random(4)
for st,name in [(base,'gmp16_base  (exactly where y was derived)'),
                (swap(base,2081,4287),'sibling: x2081 -> x4287'),
                (swap(base,2081,13195),'sibling: x2081 -> x13195')]:
    forwardp_frozen(st,set())
    bd={a:evalp(L.polys[a],st) for a in rows}
    ann=[0]*len(CERT); tot=0
    for u in rnd.sample(live,min(150,len(live))):
        v=list(st); v[u]=(v[u]+1)%P; forwardp_frozen(v,set())
        d={a:(evalp(L.polys[a],v)-bd[a])%P for a in rows}
        if not any(d.values()): continue
        tot+=1
        for j,c in enumerate(CERT):
            if sum(y*d.get(a,0) for a,y in c.items())%P==0: ann[j]+=1
    iv=[sum(y*bd[a] for a,y in c.items())%P for c in CERT]
    print(f"{name}")
    print(f"   knobs with an effect: {tot};  annihilated by each certificate: {ann}")
    print(f"   invariant values zero? {[x==0 for x in iv]}")
