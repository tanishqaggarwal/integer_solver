import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
A=load_atoms()
var_atoms=defaultdict(list)
for ai,poly in enumerate(A):
    for v in atom_vars(poly): var_atoms[v].append(ai)
def deg(poly): return max(len(m) for m in poly)
def show(ai,ml=140):
    poly=A[ai]
    terms=[f"{c}*{'*'.join(f'x_{v}' for v in m)}" if m else f"{c}" for m,c in sorted(poly.items(),key=lambda x:len(x[0]))]
    return f"a{ai}[deg{deg(poly)},{len(poly)}t]: "+(" + ".join(terms))[:ml]
for name,v in [('x_14853',14853),('x_12186',12186),('x_16742',16742)]:
    print(f"=== atoms with {name} ({len(var_atoms[v])}) ===")
    for ai in var_atoms[v]:
        print("  ",show(ai))
