import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
var_in_atom=defaultdict(list)
for ai,a in enumerate(atoms):
    vs=set()
    for mono,c in a['poly']:
        for v in mono: vs.add(v)
    for v in vs: var_in_atom[v].append(ai)

def status(v):
    isfree=v in H.freeinp
    rhs=''
    if v in H.definer:
        gi=H.definer[v]; t,r,vids=H.gates[gi]; rhs=r[:50]
    return ('FREE' if isfree else 'gate'), len(var_in_atom.get(v,[])), rhs

lvl1=[2964,26756,579,19569,24548,25442,7927,11052,
      34661,16900,5015,27289,23707,22511,  # 44342 verifier vars
      25859,25539,18312,8557,19022]         # 45677 verifier vars
for v in lvl1:
    st,na,rhs=status(v)
    print(f"x_{v}: {st} #atoms={na} rhs={rhs}")
