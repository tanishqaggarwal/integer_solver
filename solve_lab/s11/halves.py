import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, leaves
# half-tree nodes:  a=OR(5585,15347) b=OR(10083,21887) c=OR(11346,36945) d=OR(28806,30365)
HALVES=[5585,15347,10083,21887,11346,36945,28806,30365]
base=[0]*L.NVARS; fw.forward(base)
out={}
for h in HALVES:
    fr,_=leaves.cone_free(h)
    hits=[]
    for b in sorted(fr):
        v=list(base); v[b]=1; fw.forward(v)
        if v[h]==1: hits.append(b)
    out[h]=hits
    print(f"x{h}: {len(hits)} bits, first={hits[:6]}")
json.dump(out, open('halves.json','w'))
