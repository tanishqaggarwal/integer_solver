import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
import sy_build as B
p=H.p
atoms=[]; reprs={}; ateqs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); reprs[i]=d.get('repr',''); ateqs[i]=set(d.get('eqs',[]))
def ev(poly,v):
    s=0
    for m,c in poly:
        t=c
        for var in m: t*=v[var]
        s+=t
    return s
# Test1: does PURE agentA break the 16 when x_4432 moves by p?
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward(); F_agentA=set(H.fails())
H.val[9413]=1; H.val[4432]=vA.get(4432,0)+p; H.forward()
broke_agentA=set(H.fails())-F_agentA
print('PURE agentA (x_4287=0): moving x_4432 by p breaks:', sorted(broke_agentA))
# Test2: regime11 broken eqs -> their nonzero atoms
B.regime11()
base=H.val[:]
H.val[9413]=1; H.val[4432]=base[4432]+p; H.forward()
sixteen=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
# find nonzero atoms that live in these eqs
nzatoms=set()
for ai,poly in enumerate(atoms):
    if ateqs[ai] & set(sixteen):
        if ev(poly,H.val)!=0: nzatoms.add(ai)
print('nonzero atoms touching the 16 eqs:',len(nzatoms))
for ai in sorted(nzatoms):
    print(f'  atom {ai} eqs={sorted(ateqs[ai]&set(sixteen))} :: {reprs[ai][:75]}')
