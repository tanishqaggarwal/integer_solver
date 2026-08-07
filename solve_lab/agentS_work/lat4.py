"""Run the exact lattice-membership test at every distinct boolean image point."""
import sys, json, collections, pickle, time, math
sys.path.insert(0,'.')
import common as C, lat2, lat3, lattice as L
import harness as H, engine as E, fast, intsolve
P=C.P
ex=pickle.load(open('image.pkl','rb'))['ex']
print("distinct boolean image points:",len(ex))
res=[]
for i,(key,a) in enumerate(sorted(ex.items(), key=lambda kv:-len(kv[1]))):
    seed=dict(C.BASE); seed.update({int(k):int(v) for k,v in a.items()})
    try:
        r=lat3.analyse(seed,'img%d(%d bits set)'%(i,sum(1 for v in a.values() if v)))
    except Exception as e:
        print("  ERR",type(e).__name__,e); continue
    res.append((key,r is not None and r[0] is not None))
print("\nany FULL SOLVE?", any(x[1] for x in res))
