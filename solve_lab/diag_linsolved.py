import json
from collections import defaultdict
p=2**256-2**32-977
d=json.load(open('linear_solved.json'))
val=defaultdict(int)
for k,vv in d.items(): val[int(k[2:])]=int(vv)
atoms=[json.loads(l) for l in open('atoms/poly_atoms.jsonl')]
def atomval(i):
    s=0
    for vs,c in atoms[i]['poly']:
        m=c
        for x in vs: m*=val[x]
        s+=m
    return s
nz=[(i,atomval(i)) for i in range(len(atoms)) if atomval(i)!=0]
print(f"nonzero atoms: {len(nz)}")
# classify by degree and %p
bykind=defaultdict(int)
for i,v in nz:
    deg=max(len(vs) for vs,_ in atoms[i]['poly'])
    bykind[(deg, v%p==0)]+=1
for k,c in sorted(bykind.items()): print(f"  deg{k[0]}, {'p-mult' if k[1] else 'sub-p'}: {c}")
# are the nonzero deg-2 atoms DEFINITIONS (x_out - x_a*x_b) that forward should satisfy? check a few
defs=0; cons=0
for i,v in nz:
    a=atoms[i]; poly=a['poly']
    # definition: one lone var term = product term (x_out - x_a*x_b)
    lone=[vs[0] for vs,c in poly if len(vs)==1 and abs(c)==1]
    prod=[vs for vs,c in poly if len(vs)==2]
    if len(poly)==2 and len(lone)==1 and len(prod)==1: defs+=1
    else: cons+=1
print(f"of nonzero: {defs} look like pure definitions, {cons} are constraints/mixed")
print("sample nonzero atoms:")
for i,v in nz[:12]:
    print(f"  atom{i} (%p={'0' if v%p==0 else 'nz'}): {atoms[i]['repr'][:70]}")
