"""Cone-restricted evaluation of invariant 5, then the full 2^18 enumeration.

inv5 depends on 8 atoms, so only their transitive input cone needs evaluating -- not all 31,475
gates.  If the cone is small this turns a 6-hour enumeration into minutes.
"""
import sys, os, json, time, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
import fw
from gmp1 import evalp, forwardp, solvep
from bits5 import msg
P=L.P; sys.set_int_max_str_digits(400000)
CERT=[{int(k):int(v) for k,v in c.items()} for c in json.load(open(os.path.join(HERE,'data','certs.json')))]
C5=CERT[5]
# transitive input cone of the certificate's atoms
need=set()
frontier=set()
for a in C5: frontier |= set(L.avars[a])
while frontier:
    u=frontier.pop()
    if u in need: continue
    need.add(u)
    d=L.definer.get(u)
    if d is not None: frontier |= (set(L.avars[d])-need)
print(f"inv5 cone: {len(need)} variables of {L.NVARS}")
ORDER=[c for c in fw.ORDER if any(u in need for u in c)]
print(f"SCC components to evaluate: {len(ORDER)} of {len(fw.ORDER)}")
def fp_cone(v):
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
    return v
def inv5(v): return sum(y*evalp(L.polys[a],v) for a,y in C5.items())%P
ACT=[1530,2081,4287,6821,11368,12054,13195,14808,16586,17406,17760,21074,22351,22562,23751,24365,34974,35979]
# validate the cone evaluator against the full one
tpl=msg({24601})
for S in [(),(2081,),(4287,13195),(1530,17406,22562)]:
    full=msg({24601}|set(S)); a=inv5(full)
    v=list(tpl)
    D=json.load(open(os.path.join(HERE,'data','bits1.json')))
    PIN={int(k):[(x[0],x[1],int(x[2])) for x in val] for k,val in D['pins'].items()}
    for b in ACT:
        v[b]=1 if b in S else 0
        for _a,x,C in PIN[b]: v[x]=(C%P) if b in S else 0
    fp_cone(v); b2=inv5(v)
    print(f"   subset {S}: full={str(a)[:14]}.. cone={str(b2)[:14]}.. match={a==b2}")
t0=time.time()
for _ in range(200): fp_cone(list(tpl))
print(f"cone eval: {(time.time()-t0)/200*1000:.1f} ms  -> 2^18 in {(time.time()-t0)/200*262144/60:.0f} min")
