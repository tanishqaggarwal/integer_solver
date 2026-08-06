#!/usr/bin/env python3
"""Get Q = sqrt(atom 45276) and test whether x_19964=CONST1, x_2099=CONST2 is admissible in Q=0
by activating x_31861 + private slacks, keeping other vars at agentA."""
import heal_harness as H
from propagate import load_atoms, atom_vars
p=H.p
atoms=load_atoms()
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
# try_sqrt: recover degree-2 root of the degree-4 square atom
try:
    from check_square import try_sqrt
    Q=try_sqrt(atoms[45276])
    print("got Q, terms:",len(Q) if Q else None)
except Exception as e:
    print("try_sqrt failed:",e); Q=None
# Evaluate Q at agentA (should be 0)
def ev(poly,v):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
if Q:
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    H.forward()
    print("Q at agentA =", ev(Q,H.val)%p)
    # dQ/dx_31861 (coeff of x_31861 in Q, i.e. sum of its partner values)
    c31861=0
    for m,c in Q.items():
        if 31861 in m:
            if len(m)==1: c31861+=c
            elif len(m)==2:
                o=m[0] if m[1]==31861 else m[1]
                c31861+=c*H.val[o]
    print("dQ/dx_31861 at agentA =", c31861%p, "(nonzero => x_31861 can move Q)")
    # What does Q pin? coeff of x_19964, x_2099 in Q (linear sensitivity)
    for tv in [19964,2099,6418,12553]:
        cc=0
        for m,c in Q.items():
            if tv in m:
                if len(m)==1: cc+=c
                elif len(m)==2:
                    o=m[0] if m[1]==tv else m[1]; cc+=c*H.val[o]
        print(f"  dQ/dx_{tv} = {cc%p}")
