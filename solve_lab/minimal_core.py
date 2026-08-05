import heal_harness as H
from collections import defaultdict
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
F0=set(H.fails())
# is x_16742 in x_24908's cone? (circularity)
gdef={}
import json
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
seen=set(); st=[24908]
while st:
    v=st.pop()
    if v in seen: continue
    seen.add(v)
    if v in gdef:
        for u in gdef[v][1]: st.append(u)
print(f"x_16742 in x_24908 cone (circular)? {16742 in seen}")
print(f"x_14853 in x_24908 cone? {14853 in seen}, x_12186 in cone? {12186 in seen}")
LOADS=[11150,25739,37758]
print(f"BEFORE: L mod p = {[V[n]%p for n in LOADS]}")
print(f"  x_29322%p={(V[14853]-V[12186])%p}, x_3558%p={(V[24908]-V[16742])%p}")
# minimal core fix
V[14853]=V[12186]  # x_29322=0 exactly
V[16742]=(V[16742]//p)*p + V[24908]%p  # x_16742 residue = x_24908 residue
H.forward()
print(f"AFTER minimal fix: x_29322%p={(V[14853]-V[12186])%p}, x_3558%p={(V[24908]-V[16742])%p}")
print(f"  L mod p = {[V[n]%p for n in LOADS]}")
print(f"  loads all zero mod p? {all(V[n]%p==0 for n in LOADS)}")
F=set(H.fails())
print(f"  fails: {len(F)} (was {len(F0)}); newly broken: {len(F-F0)}")
# which value*value gates broke? find nonzero atoms
atoms=[json.loads(l) for l in open('atoms/poly_atoms.jsonl')]
def av(i):
    s=0
    for vs,c in atoms[i]['poly']:
        m=c
        for x in vs: m*=V[x]
        s+=m
    return s
nz=[i for i in range(len(atoms)) if av(i)!=0]
print(f"  nonzero atoms: {len(nz)}")
