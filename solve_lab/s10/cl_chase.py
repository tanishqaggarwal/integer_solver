"""CL: from the K1 gauge state (score 39009, nz [29539,31672,40826]), chase a31672."""
import os, sys, json, collections, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)
K1FREE = [8778, 14623, 16742, 31339, 33462]
d1 = (v0[24548] - v0[14623]) % P

v1 = list(v0)
for u in K1FREE: v1[u] = v1[u] + d1
ad.fwd(v1, rounds=10)
for a in [a for a in range(L.NA) if L.evalpoly(L.polys[a], v1) and a not in atom_out]:
    E.absorb(v1, a)
av1, nz1, S1, bad1 = E.stats(v1)
print(f'K1 gauge state: score {S1} nz {nz1}')
T.save(v1, os.path.join(HERE,'cl_gauge1.json'))

print('\n=== a31672 anatomy ===')
print(L.atom_src[31672])
for w in sorted(L.avars[31672]):
    print(f'   x_{w:<6} {"FREE" if w in FREE else "a"+str(definer[w]):<8} val={str(v1[w])[:34]} consumers={len(L.var_atoms[w])}')
print(f'x_24601: {"FREE" if 24601 in FREE else "gate a"+str(definer.get(24601))}  val={v1[24601]}')
for a in sorted(L.var_atoms[24601]):
    print(f'   consumer a{a:<6} {"CHECK" if a not in atom_out else "gate->x_"+str(atom_out[a][1]):<16} eqs={len(L.atom2eq.get(a,{}))} {L.atom_src[a][:110]}')

# repairs of a31672
vm1 = [x % P for x in v1]
g = ad.grad(31672, vm1)
print(f'\na31672 grad support {len(g)}: {sorted(g)}')
rows=[]
for u in sorted(g):
    if u in E.FORBID: continue
    w = E.newton(v1, 31672, u, vm=vm1, g=g, av=av1)
    if w is None: continue
    av, nz, s, bad = E.stats(w)
    rows.append((s,u,sorted(nz)))
rows.sort(key=lambda t:-t[0])
for s,u,nz in rows[:15]:
    print(f'  newton x_{u:<6} -> score {s} ({s-S1:+d})  nz={nz[:10]}')

# selector kill: set x_24601 = 0
for sel in (24601,):
    if sel in FREE:
        w = list(v1); w[sel] = 0
        ad.fwd(w, rounds=10)
        for a in [a for a in range(L.NA) if L.evalpoly(L.polys[a], w) and a not in atom_out]:
            E.absorb(w, a)
        av, nz, s, bad = E.stats(w)
        print(f'\nset x_{sel}=0 -> score {s} ({s-S1:+d}) nz={nz[:14]}')
        T.save(w, os.path.join(HERE,'cl_sel0.json'))
