import json
from collections import defaultdict
from propagate import load_atoms, atom_vars
A=load_atoms()
def show(ai):
    poly=A[ai]
    terms=[]
    for m,c in poly.items():
        if m==(): terms.append(f"{c}")
        else:
            mono="*".join(f"x_{v}" for v in m)
            terms.append(f"{c}*{mono}")
    return f"atom {ai}: " + " + ".join(terms)
for ai in [20862,20864,42669,44342,45677]:
    print(show(ai)[:400])
    print("   vars:", sorted(atom_vars(A[ai])))
