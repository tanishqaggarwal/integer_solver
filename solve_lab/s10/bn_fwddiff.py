import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
BEST = os.path.join(LAB,'best','new_instance_partial_39026.json')
v = L.load(BEST)
w = list(v); ad.fwd(w,1)
diff = [u for u in range(L.NVARS) if v[u]!=w[u]]
print('vars changed by fwd:', len(diff), diff[:40])
FREESET=set(ad.FREE)
print('  of which free:', sum(1 for u in diff if u in FREESET))
# which gate atoms are nonzero in v (i.e. broken gates)?
av=L.all_atom_values(v)
nzg=[a for a in range(L.NA) if av[a] and a in L.atom_out]
print('nonzero GATE atoms in best:', nzg)
print('nonzero CHECK atoms in best:', [a for a in range(L.NA) if av[a] and a not in L.atom_out])
# why does fwd change vars then?
for u in diff[:20]:
    d=L.definer.get(u)
    print(f'  x_{u}: v={str(v[u])[:30]} w={str(w[u])[:30]} definer=a{d} atomval={av[d] if d is not None else None}')
# repeat fwd idempotence from w
w2=list(w); ad.fwd(w2,1)
print('fwd idempotent from w:', sum(1 for u in range(L.NVARS) if w[u]!=w2[u]))
av2=L.all_atom_values(w); print('score(w)', L.NEQ-len(L.failing_eqs(av2)),
      'nonzero atoms', sum(1 for a in range(L.NA) if av2[a]))
