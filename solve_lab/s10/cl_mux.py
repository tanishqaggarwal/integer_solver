"""CL step 2b: mux structure at the top of each cone + consumers of the RHS free inputs."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v = L.load(os.path.join(HERE,'mod9118_0.json'))
av = L.all_atom_values(v)

def show(u, depth=0, seen=None, maxd=4):
    if seen is None: seen=set()
    pad='  '*depth
    a = definer.get(u)
    tag = 'FREE' if u in FREE else f'a{a}'
    val = v[u]
    vs = str(val) if abs(val)<10**12 else f'{str(val)[:14]}..({len(str(abs(val)))}d)'
    src = L.atom_src[a][:110] if a is not None else ''
    print(f'{pad}x_{u} [{tag}] = {vs}   {src}')
    if depth>=maxd or a is None or u in seen: return
    seen.add(u)
    for w in sorted(L.avars[a]):
        if w!=u: show(w, depth+1, seen, maxd)

for top in (27522, 1308):
    print(f'\n########## TOP OF CONE x_{top} ##########')
    show(top, maxd=3)

print('\n########## KEY VARIABLES ##########')
for u in (14623, 27522, 36864, 14853, 1308, 29967):
    a = definer.get(u)
    print(f'x_{u}: {"FREE" if u in FREE else f"def by a{a}: "+L.atom_src[a][:120]}')
    print(f'    value = {v[u]}   (bits {v[u].bit_length()})  mod p = {v[u]%P}')
    print(f'    consumers (atoms mentioning it): {len(L.var_atoms[u])}')

print('\n########## CONSUMERS OF x_14623 / x_14853 ##########')
for u in (14623, 14853):
    print(f'\n--- x_{u} (value {v[u]}, {v[u].bit_length()} bits) ---')
    for a in sorted(L.var_atoms[u]):
        o = atom_out.get(a)
        eqs = sorted(L.atom2eq.get(a,{}))
        print(f'  a{a:<6} out={o if o is None else "x_"+str(o[1])}  val={av[a]}  eqs={len(eqs)}{eqs if len(eqs)<14 else ""}')
        print(f'         {L.atom_src[a][:170]}')
