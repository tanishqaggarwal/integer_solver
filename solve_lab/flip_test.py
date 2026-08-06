import heal_harness as H
p=H.p
ATOMS=[]
import json
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: ATOMS.append([(tuple(m),c) for m,c in json.loads(line)['poly']])
def av(i):
    s=0
    for m,c in ATOMS[i]:
        tt=c%p
        for v in m: tt=tt*H.val[v]%p
        s=(s+tt)%p
    return s
def nzcount():
    return sum(1 for i in range(len(ATOMS)) if av(i)!=0)
for start in ['new_instance_partial_39022','new_instance_partial_39021']:
    d=H.loadd('best/'+start+'.json')
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    base=(len(H.fails()),nzcount())
    print(f"\n{start}: fails={base[0]} nonzero_atoms={base[1]}  x_2081={H.val[2081]} x_24601={H.val[24601]}")
    for bit in [2081,24601]:
        for newval in [0,1]:
            for v in H.freeinp: H.val[v]=d.get(v,0)
            H.val[bit]=newval; H.forward()
            print(f"   set x_{bit}={newval}: fails={len(H.fails())} nonzero_atoms={nzcount()}")
