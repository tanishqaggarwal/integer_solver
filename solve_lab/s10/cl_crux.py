"""CL: the exact constraint web of cluster 1 = {a7930, a21617, a29539, a33796}."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)

CL1 = [7930, 21617, 29539, 33796]
CL2 = [2423, 26731, 33929]
print('===== cluster atoms =====')
for a in CL1+CL2:
    print(f'a{a:<6} val={"0" if av0[a]==0 else "NONZERO"} eqs={len(L.atom2eq.get(a,{}))}  {L.atom_src[a][:130]}')

def cone(u):
    seen=set(); st=[u]
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w)
        a=definer.get(w)
        if a is None: continue
        for z in L.avars[a]:
            if z!=w and z not in seen: st.append(z)
    return seen

print('\n===== cone(x_25442) : the other side of a7930 =====')
C = cone(25442); F = sorted(u for u in C if u in FREE)
print(f'{len(C)} vars, {len(F)} free: ' + ', '.join(f'x_{u}({v0[u].bit_length()}b)' for u in F))
order = [t for t in ad.ORDER if t in C and t in definer]
for t in order:
    a = definer[t]
    print(f'   x_{t:<6} <- a{a:<6}: {L.atom_src[a][:120]}   = {str(v0[t])[:26]}')
print(f'x_25442 = {v0[25442]}  ({v0[25442].bit_length()} bits)  mod p = {v0[25442]%P}')
print(f'x_24548 mod p = {v0[24548]%P}')
print(f'K1 = x_14623 mod p = {v0[14623]%P}')

vm=[x%P for x in v0]
for a in CL1+CL2:
    g = ad.grad(a, vm)
    print(f'\na{a}: grad support {len(g)} : {sorted(g)[:30]}')
