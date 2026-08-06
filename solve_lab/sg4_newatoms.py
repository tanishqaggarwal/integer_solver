import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
A=load_atoms()
mem=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): mem[v]+=1
def show(ai,ml=220):
    poly=A[ai]
    deg=max(len(m) for m in poly)
    terms=[f"{c}*{'*'.join(f'x_{v}' for v in m)}" if m else f"{c}" for m,c in sorted(poly.items(),key=lambda x:(len(x[0]),x[0]))]
    return f"a{ai}[deg{deg},{len(poly)}t]: "+(" + ".join(terms))[:ml]
# the new nonzero atoms
for ai in [3269,3271,17897,20866,20868,34232,36874,40492,41740,45603]:
    print(show(ai))
    print("   vars:",sorted(atom_vars(A[ai])))
