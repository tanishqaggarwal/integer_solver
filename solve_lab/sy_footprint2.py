import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
atoms=[]; ateqs={}; reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); ateqs[i]=set(d.get('eqs',[])); reprs[i]=d.get('repr','')
# which atoms contain a given var
from collections import defaultdict
var2atoms=defaultdict(list)
for ai,poly in enumerate(atoms):
    vs=set()
    for m,c in poly: vs.update(m)
    for v in vs: var2atoms[v].append(ai)
for kn in [24548,11052,7927,25442]:
    ats=var2atoms[kn]
    tag='FREE' if kn in H.freeinp else 'gate'
    print(f'x_{kn} ({tag}): in {len(ats)} atoms: {ats[:8]}')
    for ai in ats[:6]:
        print(f'    atom {ai}: {reprs[ai][:70]}')
