import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from newton_modp import score, FORBID
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
def pr(a, n=200):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'
print('a3575:', pr(a := 3575))
for w in sorted(set(L.avars[3575])):
    print(f'   x_{w:<7} {"FREE" if w in FREE else "gate a"+str(definer[w]):<12} '
          f'consumers {len(L.var_atoms[w])}')
v = L.load(os.path.join(HERE, 'crt_39009.json'))
av = L.all_atom_values(v)
print(f'\nstate score {score(v)}  nonzero {[b for b in range(L.NA) if av[b]]}')
tgt = T.solve_lin(3576, 26777, v)
print(f'x_26777 target {str(tgt)[:34]}  (current {str(v[26777])[:34]})')
vv = list(v); vv[26777] = tgt
for u in sorted(set(L.avars[3575])):
    if u == 26777 or u in FORBID: continue
    nv = T.solve_lin(3575, u, vv)
    if nv is None:
        print(f'   via x_{u}: not solvable'); continue
    cands = [(u, nv)] if u in FREE else []
    if u not in FREE:
        d = definer.get(u)
        if d is not None:
            vvv = list(v); vvv[u] = nv
            for z in sorted(set(L.avars[d])):
                if z == u or z not in FREE or z in FORBID: continue
                zz = T.solve_lin(d, z, vvv)
                if zz is not None: cands.append((z, zz))
    for z, zv in cands:
        tr = list(v); tr[z] = zv
        ad.fwd(tr, rounds=6)
        at = L.all_atom_values(tr)
        s = score(tr)
        print(f'   via x_{u} -> set x_{z}: score {s}  '
              f'nonzero {[b for b in range(L.NA) if at[b]]}')
        if s > score(v):
            T.save(tr, os.path.join(HERE, f'crt2_{s}.json'))
            print(f'      *** saved crt2_{s}.json')
