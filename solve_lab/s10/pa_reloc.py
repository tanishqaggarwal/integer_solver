"""Variant B: let the canonical frame REPAIR the delivered witness's broken gates and
see where the residual relocates to.  Then price alternative relocations."""
import os, sys, collections, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P=ad.P
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))

def report(w,tag):
    av=L.all_atom_values(w)
    nz=[a for a in range(L.NA) if av[a]]
    E=set()
    for a in nz: E|=set(L.atom2eq[a])
    f=L.failing_eqs(av)
    print(f'[{tag}] |S|={len(nz)} |E|={len(E)} failing={len(f)} score={L.NEQ-len(f)}')
    print(f'    S={nz[:30]}')
    print(f'    fail={f[:25]}')
    return av,nz,f

report(v0,'delivered')
w=list(v0); ad.fwd(w,4)
report(w,'canonical fwd from delivered')

# repair each broken gate individually (ripple, canonical frame elsewhere)
BROKEN={22229:7068,22230:28730,35758:29854,35761:31864,35762:642}
for a,t in BROKEN.items():
    w=list(v0)
    nv=T.solve_lin(a,t,w)
    if nv is None:
        print(f'repair a{a} via x_{t}: no integer solution'); continue
    blk=set(BROKEN)-{a}
    L.ripple(w,{t:nv},block=blk)
    report(w,f'repair a{a} (x_{t}={str(nv)[:14]}..)')
