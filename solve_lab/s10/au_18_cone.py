import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def cone(u):
    seen=set(); st=[u]
    while st:
        t=st.pop()
        if t in seen: continue
        seen.add(t)
        a=L.definer.get(t)
        if a is None: continue
        for x in L.avars[a]:
            if x!=t and x not in seen: st.append(x)
    return seen
for target in (2099, 19964):
    c=cone(target)
    print(f'===== cone of x_{target}: {len(c)} vars =====')
    # topological print
    order=[t for t in L.topo if t in c]
    order = [t for t in c if t not in set(order)] + order
    for t in order:
        a=L.definer.get(t)
        if a is None:
            print(f'  FREE x_{t} = {str(v[t])[:50]}  (nvars_using={len(L.var_atoms[t])})')
        else:
            print(f'  x_{t} <- a{a}: {L.atom_src[a][:110]}   val={str(v[t])[:34]}')
    print()
