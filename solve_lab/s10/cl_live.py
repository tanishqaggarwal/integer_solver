"""CL: collapse the live path of x_27522 / x_1308 and identify the pinning constraints."""
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

def live_chain(u, depth=0, seen=None, out=None):
    """follow the variables that actually carry the value (nonzero contribution)."""
    if seen is None: seen=set()
    if out is None: out=[]
    if u in seen: return out
    seen.add(u)
    a = definer.get(u)
    pad = '  '*depth
    if a is None:
        out.append((depth, u, 'FREE', v0[u]))
        print(f'{pad}x_{u} = FREE  ({v0[u].bit_length()} bits)  consumers={len(L.var_atoms[u])}')
        return out
    # which monomials of the definer are nonzero?
    live = []
    for m, c in L.polys[a].items():
        if u in m: continue
        t = c
        for z in m: t *= v0[z]
        if t: live.append((m, c, t))
    print(f'{pad}x_{u} <- a{a}: {L.atom_src[a][:100]}   [{len(live)} live monomials]')
    for m, c, t in live:
        for z in m:
            if v0[z] != 0:
                live_chain(z, depth+1, seen, out)
    return out

for top in (27522, 1308):
    print(f'\n########## LIVE CHAIN of x_{top} (value {v0[top].bit_length()} bits) ##########')
    live_chain(top)

print('\n########## the two "equal free inputs" constraints ##########')
for (chk, lhs, rhs) in [(21617, 24548, 14623), (29539, None, 14853)]:
    print(f'\na{chk}: {L.atom_src[chk][:120]}')
    if lhs: print(f'   x_{lhs} mod p = {v0[lhs]%P}')
    print(f'   x_{rhs} mod p = {v0[rhs]%P}')

for u in (24548, 14623, 14853):
    print(f'\n===== consumers of x_{u} (free, {v0[u].bit_length()} bits, val mod p = {v0[u]%P}) =====')
    for a in sorted(L.var_atoms[u]):
        o = atom_out.get(a)
        eqs = sorted(L.atom2eq.get(a,{}))
        print(f'  a{a:<6} {"CHECK" if o is None else "gate->x_"+str(o[1]):<14} val={"0" if av0[a]==0 else str(av0[a])[:18]+"..."} eqs={len(eqs)}')
        print(f'     {L.atom_src[a][:190]}')
