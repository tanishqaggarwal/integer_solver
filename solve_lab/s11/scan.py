import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, leaves
base=[0]*L.NVARS
fw.forward(base)
print("base U,V =", base[7715], base[34554])
res={}
for tgt,name in [(8599,'a'),(21839,'b'),(7304,'c'),(25956,'d')]:
    fr,_=leaves.cone_free(tgt)
    hits=[]
    for b in sorted(fr):
        v=list(base); v[b]=1
        fw.forward(v)
        if v[tgt]==1: hits.append(b)
    res[name]=hits
    print(f"x{tgt} ({name}): {len(hits)}/{len(fr)} single bits set it -> {hits[:12]}")
import json; json.dump(res, open('bitsets.json','w'))
