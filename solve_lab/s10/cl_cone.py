"""CL step 2: trace the ancestor cones of x_27522 and x_1308, classify each gate."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v]

def cone(u):
    seen=set(); st=[u]
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w)
        a = definer.get(w)
        if a is None: continue
        for z in L.avars[a]:
            if z!=w and z not in seen: st.append(z)
    return seen

def classify(a, t, v):
    """classify the gate atom a defining var t."""
    Pp = L.polys[a]
    # coefficient structure: a = c_t * x_t + rest
    r = T.lin_parts(a, t, v)
    terms=[]
    for m,c in Pp.items():
        terms.append((m,c))
    deg = max(len(m) for m in Pp)
    nterm = len(Pp)
    # is any monomial containing t of degree >1?  (t multiplied by a var)
    tmul = [m for m in Pp if t in m and len(m)>1]
    return deg, nterm, tmul

def dump(u, name, maxn=100000):
    C = cone(u)
    free = sorted(x for x in C if x in FREE)
    comp = [x for x in C if x not in FREE]
    print(f'\n===== cone({name}=x_{u}): {len(C)} vars, {len(free)} free inputs, {len(comp)} computed =====')
    # topo order restricted
    order = [t for t in ad.ORDER if t in C and t in definer]
    kinds = collections.Counter()
    for t in order:
        a = definer[t]
        Pp = L.polys[a]
        deg = max(len(m) for m in Pp)
        nterm = len(Pp)
        vs = sorted(L.avars[a])
        # monomials mentioning t
        tm = [(m,c) for m,c in Pp.items() if t in m]
        others = [(m,c) for m,c in Pp.items() if t not in m]
        if len(tm)==1 and len(tm[0][0])==1:
            # a = c*x_t + rest -> x_t = -rest/c
            if all(len(m)==1 for m,_ in others):
                k = 'LIN'         # pure linear combination
            elif all(len(m)<=2 for m,_ in others):
                k = 'BILIN'
            else:
                k = 'POLY'
        else:
            k = 'TMUL'
        kinds[k]+=1
    print('  gate kinds:', dict(kinds))
    return C, free, order

C1, F1, O1 = dump(27522, 'x_27522')
C2, F2, O2 = dump(1308, 'x_1308')
print(f'\ncone overlap vars: {len(C1&C2)}   free overlap: {len(set(F1)&set(F2))}')

# print the definition chain of the SMALL cone in full
print('\n===== FULL CHAIN of cone(x_27522) =====')
for t in O1:
    a = definer[t]
    src = L.atom_src[a]
    print(f'  x_{t:<6} <- a{a:<6} : {src[:160]}')
