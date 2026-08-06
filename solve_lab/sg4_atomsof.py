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
def evalp(poly):
    s=0
    for m,c in poly.items():
        t=c
        for v in m: t*=V[v]
        s+=t
    return s
def show(ai,maxlen=180):
    poly=A[ai]
    deg=max(len(m) for m in poly)
    terms=[]
    for m,c in sorted(poly.items(),key=lambda x:(len(x[0]),x[0])):
        if m==(): terms.append(f"{c}")
        else: terms.append(f"{c}*{'*'.join(f'x_{v}' for v in m)}")
    s=" + ".join(terms)
    return f"a{ai}[deg{deg},{len(poly)}t,r={evalp(poly)}]: {s[:maxlen]}"
for name,v in [('x_4432',4432),('x_7068',7068)]:
    print(f"===== atoms containing {name} ({len(var_atoms[v])}) =====")
    for ai in var_atoms[v]:
        print(" ",show(ai))
