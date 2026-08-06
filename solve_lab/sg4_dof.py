import heal_harness as H
import json
from collections import defaultdict
p=H.p
# Which core vars are free inputs?
core_vars=[14853,12186,24908,16742,29322,3558,4432,7068,19964,2099,642,28730,17325,9413,28599,17499]
print("var    free?   #eqs   #atoms")
# count eq appearances
eqcount=defaultdict(int)
for vs in H.eqvars:
    for v in vs: eqcount[v]+=1
from propagate import load_atoms, atom_vars
A=load_atoms()
atomcount=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): atomcount[v]+=1
for v in core_vars:
    print(f"x_{v:<7} {str(v in H.freeinp):<6} {eqcount[v]:<6} {atomcount[v]}")
print("---ancestors (free inputs each core var depends on)---")
for v in [29322,3558,14853,12186,24908,16742]:
    a=H.anc.get(v,set())
    print(f"x_{v}: {len(a)} free ancestors:", sorted(a)[:20], '...' if len(a)>20 else '')
