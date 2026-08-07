"""Honest scope of the equation-level theorem: at each window level, how many variables
appear in the window's atoms, how many are knobs, and what the excluded ones touch."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import pick_knobs
P=env.P
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fe=L.failing_eqs(L.all_atom_values(v))
A=set(a for e in fe for a in L.eq_atoms[e][2])
for lev in range(7):
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    V=set(u for a in A for u in L.avars[a])
    K=set(pick_knobs(v,A))
    excl=V-K
    # why excluded: touches an atom outside A, or dropped by the linearity filter
    out_atom=[u for u in excl if any(x not in A for x in L.var_atoms[u])]
    lin_drop=[u for u in excl if all(x in A for x in L.var_atoms[u])]
    print('lev%d atoms=%-5d eqs=%-5d vars=%-5d knobs=%-4d excluded=%-4d (touch an atom outside the window: %d ; dropped to keep linearity: %d)'%(
        lev,len(A),len(R),len(V),len(K),len(excl),len(out_atom),len(lin_drop)),flush=True)
    A=set(a for e in R for a in L.eq_atoms[e][2])
