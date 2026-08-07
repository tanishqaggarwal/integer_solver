"""Image of the map (all 256 cluster booleans) -> (a20215 mod p, a28647 mod p, ...) at cfg0 base."""
import sys, json, collections, itertools, random, time, pickle
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
base=dict(C.BASE)
v0=E.forward(base); bad0=E.badatoms(v0)
CAND=C.cluster_cone(); BOOLS=[f for f in CAND if C.isbool(f)]
ROWS=C.ROWS
ON=[f for f in BOOLS if v0[f]==1]
print("bools",len(BOOLS),"currently ON:",ON)
def g(assign):
    v,_=fast.apply_delta(v0,assign)
    ns={'v':v,'__builtins__':{}}
    return tuple(eval(H.acodes[a],ns)%P for a in ROWS), v
rnd=random.Random(3)
img=collections.Counter(); ex={}
t0=time.time(); N=int(sys.argv[1]) if len(sys.argv)>1 else 400
for it in range(N):
    k=rnd.choice([1,2,3,5,8,13,21,34,55,89,128,200])
    S=rnd.sample(BOOLS,k)
    a={f:1 for f in S}
    for f in ON:
        if f not in a and rnd.random()<0.5: a[f]=0
    key,v=g(a)
    img[key]+=1
    if key not in ex: ex[key]=dict(a)
    if key[3]==0 and key[4]==0:
        print("*** HIT (a20215,a28647) both 0 mod p ***",flush=True)
        json.dump({str(x):str(y) for x,y in a.items()},open('image_hit.json','w'))
print(f"{N} samples in {time.time()-t0:.0f}s -> {len(img)} distinct 5-tuples mod p")
for key,c in img.most_common(20):
    print(f"  x{c:4d}  20215={str(key[3])[:14]}.. 28647={str(key[4])[:14]}..  7389={str(key[0])[:10]} 10187={str(key[1])[:10]} 20212={str(key[2])[:10]}")
z=collections.Counter()
for key,c in img.items():
    z[(key[3]==0,key[4]==0)]+=c
print("zero-pattern counts (20215==0, 28647==0):",dict(z))
pickle.dump({'img':dict(img),'ex':ex},open('image.pkl','wb'))
