import heal_harness as H
import json, glob
p=H.p
ATOMS=[]; reprs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line); ATOMS.append([(tuple(m),c) for m,c in dd['poly']]); reprs.append(dd.get('repr',''))
def av(i,val):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*val[v]%p
        s=(s+tt)%p
    return s
files=['new_instance_partial_39007','new_instance_partial_39013','new_instance_partial_39018',
       'new_instance_partial_39021','new_instance_partial_39022']
allnz={}
for f in files:
    d=H.loadd('best/'+f+'.json')
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    nz=[i for i in range(len(ATOMS)) if av(i,H.val)!=0]
    allnz[f]=set(nz)
    print(f"{f}: {len(H.fails())} fails, nonzero atoms({len(nz)}): {nz}")
print("\n=== union of all obstruction atoms ===")
U=set()
for s in allnz.values(): U|=s
print(sorted(U))
for a in sorted(U):
    inwhich=[f.split('_')[-1] for f in files if a in allnz[f]]
    print(f"  atom {a}: {reprs[a][:70]}  | nonzero in: {inwhich}")
