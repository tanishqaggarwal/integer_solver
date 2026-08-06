import json, re
import heal_harness as H
from collections import defaultdict
p=H.p
# 1. wire members: signed union-find over 2-term identity atoms x_a - x_b and x_a + x_b
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
# find wire class of x_26064 via identity atoms (x_a +/- x_b, no const, no quad)
par=list(range(H.NVARS)); sgn=[1]*H.NVARS
def find(x):
    if par[x]==x: return x,1
    r,s=find(par[x]); par[x]=r; sgn[x]*=s; return r,sgn[x]
def uni(a,b,s):  # x_a = s * x_b
    ra,sa=find(a); rb,sb=find(b)
    if ra==rb: return
    par[rb]=ra; sgn[rb]=sa*s*sb
for a in atoms:
    poly=a['poly']
    if all(len(vs)<=1 for vs,c in poly):  # linear
        lin=[(vs[0],c) for vs,c in poly if len(vs)==1]
        const=sum(c for vs,c in poly if len(vs)==0)
        if len(lin)==2 and const==0 and abs(lin[0][1])==1 and abs(lin[1][1])==1:
            (va,ca),(vb,cb)=lin
            uni(va,vb,-(ca*cb))  # ca*va+cb*vb=0 -> va = -(cb/ca) vb = -(ca*cb) vb since |c|=1
r26064,_=find(26064)
wire=set(v for v in range(H.NVARS) if find(v)[0]==r26064)
print(f"wire members: {len(wire)}")
# 2. extract linear-mod-p relations: for each atom, drop wire-product terms (≡0 mod p) and wire linear terms; if remaining deg<=1, it's linear mod p
# represent relation as dict var->coeff (mod p) + const
lin_rels=[]
for ai,a in enumerate(atoms):
    poly=a['poly']
    rel=defaultdict(int); const=0; ok=True
    for vs,c in poly:
        vs=tuple(vs)
        if len(vs)==0: const=(const+c)%p
        elif len(vs)==1:
            if vs[0] in wire: pass  # ≡0 mod p
            else: rel[vs[0]]=(rel[vs[0]]+c)%p
        elif len(vs)==2:
            if vs[0] in wire or vs[1] in wire: pass  # product ≡0 mod p
            else: ok=False; break  # genuine quadratic
        else: ok=False; break
    if ok and any(v%p for v in rel.values()):  # non-trivial linear relation
        lin_rels.append((ai, dict(rel), const%p))
print(f"linear-mod-p relations (non-wire): {len(lin_rels)}")
# save
json.dump({'wire':sorted(wire),'nrels':len(lin_rels)}, open('lin_meta.json','w'))
import pickle
pickle.dump(lin_rels, open('lin_rels.pkl','wb'))
print("saved lin_rels.pkl")
# quick stats: how many distinct non-wire vars appear
vs_all=set()
for ai,rel,c in lin_rels: vs_all|=set(rel)
print(f"distinct non-wire vars in linear relations: {len(vs_all)}")
