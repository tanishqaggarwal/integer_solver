#!/usr/bin/env python3
"""At wire=1, find the TRUE failing-equation obstruction (bits free). Fix G1/G2 via fine slacks.
Identify which equations fail and their dependence on the boolean bits."""
import json,pickle
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
D=pickle.load(open('wire_data.pkl','rb')); wire=D['wire']
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v=vA[:]
# set wire -> sign*1
for y,s in wire.items(): v[y]=s*1
# equations
import re
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[frozenset(int(m) for m in VAR.findall(L)) for L in lines]
ns={'v':v,'__builtins__':{}}
F=[i for i,c in enumerate(eqcode) if eval(c,ns)!=0]
print(f"wire=1 (bits at agentA, else agentA): {len(F)} eqs fail")
print("failing:",F)
# the boolean bits
BITS=[31342,32058,22473,17389,1488,28827,37384,11094,18211,875,29159,37076,14048]
# which failing eqs contain which bits / the forcing var x_26064 / G1G2 vars
gapvars={642,2099,7068,4432,19964,28730}
for i in F:
    vs=eqvars[i]
    b=sorted(vs&set(BITS))
    has26064='x_26064' if 26064 in vs else ''
    hasgap=sorted(vs&gapvars)
    print(f"  eq {i}: bits={b} {has26064} gap={hasgap}")
