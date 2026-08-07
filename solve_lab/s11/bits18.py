"""Are the certificates valid across DIFFERENT messages?

INV_y is conserved only for moves in the span of the knobs used to derive y.  y was computed at
the checkpoint's message.  If the same y still annihilates the responses at a different message,
the invariants are structural and can screen candidate messages cheaply; if not, each message
needs its own certificate and the screen is worthless.
"""
import sys, os, json, time, random
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
FREEV=[u for u in range(L.NVARS) if u not in L.definer]
rnd=random.Random(3)
sample=rnd.sample(FREEV,120)
for S,name in [({2081,24601},'checkpoint message (where y was derived)'),
               ({4287,24601},'sibling message {x4287,x24601}'),
               ({13195,24601},'sibling message {x13195,x24601}'),
               ({24601},'{x24601} (channel C)')]:
    base=msg(S)
    bd={}
    rows=set()
    for c in CERT: rows |= set(c)
    for a in rows: bd[a]=evalp(L.polys[a],base)
    ann=[0]*len(CERT); tot=0
    for u in sample:
        v=list(base); v[u]=(v[u]+1)%P; forwardp(v)
        d={a:(evalp(L.polys[a],v)-bd[a])%P for a in rows}
        if not any(d.values()): continue
        tot+=1
        for j,c in enumerate(CERT):
            if sum(y*d.get(a,0) for a,y in c.items())%P==0: ann[j]+=1
    print(f"{name}: {tot} knobs with an effect; certificate annihilates {ann} of them")
