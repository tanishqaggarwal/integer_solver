import heal_harness as H
import json
p=H.p
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
def av(i,val=None):
    val=val or H.val; s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
def atomvars(i): return set(x for m,c in ATOMS[i] for x in m)
d=H.loadd('best/new_instance_partial_39021.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F0=set(H.fails())
M1=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
M2=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
print("before:", len(F0), "fails")
print("atoms with x_22152:", [i for i in range(len(ATOMS)) if 22152 in atomvars(i)])
print("atoms with x_33462:", [i for i in range(len(ATOMS)) if 33462 in atomvars(i)])
# direct fix
H.val[22152]=M1; H.val[33462]=M2
H.forward()
F1=set(H.fails())
print("after set x_22152=M1,x_33462=M2:", len(F1), "fails; broken:",sorted(F1-F0),"fixed:",sorted(F0-F1))
nz=[i for i in range(len(ATOMS)) if av(i)!=0]
print("nonzero atoms now:", nz)
for i in nz: print(f"   atom {i}: {reprs[i][:75]}")
