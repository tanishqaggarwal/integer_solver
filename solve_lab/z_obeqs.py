import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
eq_atoms=defaultdict(list)
for ai,a in enumerate(atoms):
    for e in set(a['eqs']): eq_atoms[e].append(ai)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1; H.forward()
def ev(ai):
    s=0
    for mono,c in atoms[ai]['poly']:
        t=c
        for v in mono: t*=H.val[v]
        s+=t
    return s%p
for e in [9346,10393,14635,18629]:
    print(f"=== eq {e} ===")
    for ai in eq_atoms[e]:
        val=ev(ai)
        if val!=0:
            print(f"  atom {ai}: {atoms[ai]['repr'][:90]}  [nonzero]")
# check the pins for 13195: what are x_18623,x_6467 CONSTs and did setting help
pins=json.load(open('pinrec.json'))
p13195=[r for r in pins if r[1]==13195]
print("\npins for bit 13195:")
for r in p13195: print(f"  atom{r[0]} target=x_{r[2]} const={str(r[3])[:20]}.. coef={r[4]} handle=x_{r[5]}")
