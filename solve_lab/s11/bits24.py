"""COMPLETE enumeration of invariant 5 over all 2^18 subsets of the bits that move it.

inv5 is a function of the C-subset alone, and only 18 C bits move it, so this is exhaustive:
either 0 is attainable and we get the subset, or invariant 5 can never vanish -- which would say
something decisive about the whole instance.
"""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
import fw
from gmp1 import evalp, solvep
from bits5 import msg
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
need=set(); frontier=set()
for a in C5: frontier |= set(L.avars[a])
while frontier:
    u=frontier.pop()
    if u in need: continue
    need.add(u)
    d=L.definer.get(u)
    if d is not None: frontier |= (set(L.avars[d])-need)
ORDER=[c for c in fw.ORDER if any(u in need for u in c)]
def fp(v):
    for comp in ORDER:
        if len(comp)==1:
            u=comp[0]
            if u not in need: continue
            x=solvep(L.definer[u],u,v)
            if x is not None: v[u]=x
        else:
            for _ in range(60):
                ch=False
                for u in comp:
                    if u not in need: continue
                    x=solvep(L.definer[u],u,v)
                    if x is not None and x!=v[u]: v[u]=x; ch=True
                if not ch: break
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
ACT=[1530,2081,4287,6821,11368,12054,13195,14808,16586,17406,17760,21074,22351,22562,23751,24365,34974,35979]
D=json.load(open(os.path.join(HERE,'data','bits1.json')))
PIN={int(k):[(x[0],x[1],int(x[2])) for x in val] for k,val in D['pins'].items()}
tpl=msg({24601})
for b in ACT:
    tpl[b]=0
    for _a,x,C in PIN[b]: tpl[x]=0
N=1<<len(ACT)
t0=time.time(); vals=collections.Counter(); zeros=[]
for m in range(N):
    v=list(tpl)
    for i,b in enumerate(ACT):
        if m>>i&1:
            v[b]=1
            for _a,x,C in PIN[b]: v[x]=C%P
    fp(v)
    x=inv5(v); vals[x]+=1
    if x==0:
        S=[ACT[i] for i in range(len(ACT)) if m>>i&1]
        zeros.append(S); print(f"  *** inv5 == 0 at C-subset {S}", flush=True)
    if m%20000==0:
        print(f"   {m}/{N} distinct={len(vals)} zeros={len(zeros)} ({time.time()-t0:.0f}s)", flush=True)
print(f"COMPLETE: {N} subsets, {len(vals)} distinct inv5 values, zeros {len(zeros)} ({time.time()-t0:.0f}s)")
print("  multiplicity profile:", collections.Counter(vals.values()).most_common(6))
json.dump({'N':N,'distinct':len(vals),'zeros':zeros}, open(os.path.join(HERE,'data','bits24.json'),'w'))
