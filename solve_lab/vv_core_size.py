import json, re
import heal_harness as H
from collections import defaultdict
p=H.p
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
pinrec=json.load(open('pinrec.json'))
bits={sel for i,sel,tgt,const,coef,handle in pinrec}
# value*value gates (both factors non-wire AND non-bit)
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
VAR=re.compile(r'x_(\d+)')
vv_out=set()  # outputs of genuine value*value gates
for t,(rhs,vids) in gdef.items():
    if '*' in rhs:
        vs=[int(m) for m in VAR.findall(rhs)]
        if len(vs)==2 and vs[0] not in wire and vs[1] not in wire and vs[0] not in bits and vs[1] not in bits:
            vv_out.add(t)
print(f"genuine value*value gate outputs (non-wire non-bit factors): {len(vv_out)}")
# which value*value outputs are ancestors of the verifier-square (deg-4) atoms / core?
# backward cone from constraint atoms. constraint atoms: deg-4 squares + the 20 core eqs' atoms
# get vars in deg-4 atoms
constr_vars=set()
for a in atoms:
    deg=max(len(vs) for vs,_ in a['poly'])
    if deg>=4:
        for vs,c in a['poly']:
            for v in vs: constr_vars.add(v)
# also the core control vars
constr_vars|={14853,12186,24908,16742,3558,29322,35389,6671,11150,25739,37758}
# backward cone: all ancestors (through gate defs) of constr_vars
anc_all=set()
stack=list(constr_vars)
seen=set()
while stack:
    v=stack.pop()
    if v in seen: continue
    seen.add(v)
    if v in gdef:
        for u in gdef[v][1]:
            if u not in seen: stack.append(u)
vv_on_paths=vv_out & seen
print(f"value*value outputs in backward cone of constraints: {len(vv_on_paths)}")
# of these, how many have BOTH factors also depending on the core control vars (truly coupling)?
print(f"total vars in constraint cone: {len(seen)}")
# how many value*value gates are 'squares' x_a*x_a (like x_29356=x_29322^2)?
sq=sum(1 for t in vv_on_paths if len(set(gdef[t][1]))==1)
print(f"  of which squares (x_a^2): {sq}")
