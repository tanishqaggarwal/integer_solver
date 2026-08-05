#!/usr/bin/env python3
"""Find x_14853, x_12186 check atoms; find unchecked frees controlling x_29322 mod p with LOW ripple."""
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
import heal_harness as H
p=2**256-2**32-977
atoms=load_atoms()
CK=pickle.load(open('checked.pkl','rb')); checked=CK['checked']; checkatom=CK['checkatom']
# x_14853, x_12186 check atoms
for v in [14853,12186,16742]:
    if v in checkatom:
        ai=checkatom[v]; poly=atoms[ai]
        terms=[(m,c) for m,c in poly.items()]
        print(f"x_{v} check atom {ai}: "+" ".join(f"{c:+d}*x_{m[0]}" if m else f"{c:+d}" for m,c in terms))
# fanout (eqs touched) per free
eq_free=[]
for i in range(len(H.eqcode)):
    s=set()
    for w in H.eqvars[i]:
        if w in H.freeinp: s.add(w)
        s|=H.anc.get(w,set())
    eq_free.append(s)
free_eqs=defaultdict(int)
for i in range(len(H.eqcode)):
    for w in eq_free[i]: free_eqs[w]+=1
# x_29322 = x_14853 - x_12186; its free-ancestors (unchecked) and their fanout
anc29322=H.anc.get(29322,set())
anc3558=H.anc.get(3558,set())
print(f"\nx_29322 free-ancestors: {len(anc29322)}; x_3558: {len(anc3558)}")
unchk_29322=[(w,free_eqs[w],w in checked) for w in anc29322 if w not in checked]
unchk_29322.sort(key=lambda t:t[1])
print("x_29322 UNCHECKED ancestors (by fanout):", [(w,f) for w,f,_ in unchk_29322[:12]])
unchk_3558=[(w,free_eqs[w]) for w in anc3558 if w not in checked]
unchk_3558.sort(key=lambda t:t[1])
print("x_3558 UNCHECKED ancestors (by fanout):", unchk_3558[:12])
# lowest-fanout unchecked control for each
