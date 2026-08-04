import json, re
import heal_harness as H
from collections import defaultdict
p=H.p
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
# wire members
par2=list(range(H.NVARS))
def f3(x):
    while par2[x]!=x: par2[x]=par2[par2[x]]; x=par2[x]
    return x
atoms=[json.loads(l) for l in open('atoms/poly_atoms.jsonl')]
for a in atoms:
    poly=a['poly']
    if all(len(vs)<=1 for vs,c in poly):
        lin=[(vs[0],c) for vs,c in poly if len(vs)==1]; const=sum(c for vs,c in poly if len(vs)==0)
        if len(lin)==2 and const==0 and abs(lin[0][1])==1 and abs(lin[1][1])==1:
            ra,rb=f3(lin[0][0]),f3(lin[1][0])
            if ra!=rb: par2[rb]=ra
wire=set(v for v in range(H.NVARS) if f3(v)==f3(26064))
pinrec=json.load(open('pinrec.json')); bits={s for i,s,t,c,co,h in pinrec}
VAR=re.compile(r'x_(\d+)')
# trace x_24908 backward, mark value*value gates encountered
def trace(root, maxdepth=40):
    seen=set(); vv=[]; stack=[(root,0)]
    while stack:
        v,d=stack.pop()
        if v in seen or d>maxdepth: continue
        seen.add(v)
        if v in gdef:
            rhs,vids=gdef[v]
            if '*' in rhs:
                fs=[int(m) for m in VAR.findall(rhs)]
                if len(fs)==2 and fs[0] not in wire and fs[1] not in wire and fs[0] not in bits and fs[1] not in bits:
                    vv.append((v,rhs,fs))
            for u in vids: stack.append((u,d+1))
    return seen,vv
for name in [24908,17601,3558,29322,35389,6671]:
    seen,vv=trace(name)
    print(f"x_{name}: {len(seen)} vars in cone, {len(vv)} genuine value*value gates")
    if name==24908:
        print(f"  x_24908 def: {gdef.get(24908)}")
        print(f"  first few value*value on path: {[(t,rhs[:40]) for t,rhs,fs in vv[:5]]}")
# key: is x_24908 reachable to free inputs LINEARLY (no value*value), or through value*value?
seen,vv=trace(24908)
freeanc=seen & H.freeinp
print(f"\nx_24908: {len(freeanc)} free ancestors, {len(vv)} value*value gates on path")
print(f"  -> x_24908 is {'LINEAR in free inputs' if len(vv)==0 else 'VALUE*VALUE-coupled'}")
