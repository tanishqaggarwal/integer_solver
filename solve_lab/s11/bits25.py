"""Full sweep: do the certificates annihilate EVERY live knob at other messages?

bits19 only found 1-2 knobs touching the certificate rows out of 150 sampled, which is thin.
Sweep all 1,470 live non-bit knobs at several messages and count exactly.
"""
import sys, os, json, time, pickle, random
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
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
rows=sorted(rows)
for S,name in [({2081,24601},'checkpoint {x2081,x24601}'),
               ({4287,24601},'sibling  {x4287,x24601}'),
               ({13195,24601},'sibling  {x13195,x24601}')]:
    base=msg(S); forwardp(base)
    bd={a:evalp(L.polys[a],base) for a in rows}
    ann=[0]*len(CERT); tot=0; t0=time.time()
    for u in live:
        v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
        d={a:(evalp(L.polys[a],v)-bd[a])%P for a in rows}
        if not any(d.values()): continue
        tot+=1
        for j,c in enumerate(CERT):
            if sum(y*d.get(a,0) for a,y in c.items())%P==0: ann[j]+=1
    iv=[sum(y*bd[a] for a,y in c.items())%P for c in CERT]
    print(f"{name}: {tot} of {len(live)} live knobs touch the certificate rows; "
          f"annihilated {ann}  ({time.time()-t0:.0f}s)")
    print(f"    invariants zero? {[x==0 for x in iv]}", flush=True)
