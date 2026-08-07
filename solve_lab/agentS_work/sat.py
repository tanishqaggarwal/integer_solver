"""Rigorous saturation test: g(S) = (a20215 mod p, a28647 mod p) as a function of the bit set S."""
import sys, json, collections, itertools, random, time, pickle
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
base=dict(C.BASE)
v0=E.forward(base); bad0=E.badatoms(v0)
CAND=C.cluster_cone(); BOOLS=[f for f in CAND if C.isbool(f)]
print("bools in cluster cone:",len(BOOLS))
ns0={'v':v0,'__builtins__':{}}
def g(S):
    ch={f:1 for f in S}
    v,_=fast.apply_delta(v0,ch)
    ns={'v':v,'__builtins__':{}}
    return (eval(H.acodes[20215],ns)%P, eval(H.acodes[28647],ns)%P)
G0=g([])
print("g(empty) =",G0)
# per-bit deltas
d1={}
t0=time.time()
for f in BOOLS:
    if v0[f]==1: continue
    G=g([f]); d=((G[0]-G0[0])%P,(G[1]-G0[1])%P)
    if any(d): d1[f]=d
print("moving bits: %d  (%.0fs)"%(len(d1),time.time()-t0))
cls=collections.defaultdict(list)
for f,d in d1.items(): cls[d].append(f)
keys=sorted(cls,key=lambda k:-len(cls[k]))
print("classes:",[(len(cls[k]),str(k[0])[:10],str(k[1])[:10]) for k in keys])
rnd=random.Random(7)
print("\n--- ADDITIVITY / SATURATION TEST ---")
tests=[]
for i,k in enumerate(keys):
    m=cls[k]
    if len(m)>=2: tests.append(('within cls%d'%i, rnd.sample(m,2)))
    if len(m)>=3: tests.append(('within cls%d x3'%i, rnd.sample(m,3)))
for i in range(len(keys)):
    for j in range(i+1,len(keys)):
        tests.append(('cross cls%d+cls%d'%(i,j),[cls[keys[i]][0],cls[keys[j]][0]]))
for tag,S in tests:
    G=g(S)
    act=((G[0]-G0[0])%P,(G[1]-G0[1])%P)
    add=(sum(d1[f][0] for f in S)%P, sum(d1[f][1] for f in S)%P)
    sat=None
    # is it equal to the delta of a single member of some class present?
    singles={f:d1[f] for f in S}
    verdict='ADDITIVE' if act==add else ('SATURATED(=one member)' if act in singles.values() else 'OTHER')
    print(f"  {tag:22s} S={S}  {verdict}")
    if verdict=='OTHER':
        print(f"      act=({str(act[0])[:22]},{str(act[1])[:22]}) add=({str(add[0])[:22]},{str(add[1])[:22]})")
pickle.dump({'G0':G0,'d1':d1,'cls':dict(cls)},open('sat_cfg0.pkl','wb'))
