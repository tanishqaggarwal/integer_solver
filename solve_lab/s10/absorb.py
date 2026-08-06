"""S10 step 97: after the mod-p Newton step, make the p-handle absorb a/p over Z."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from newton_modp import newton_moves, score, FORBID
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
def pr(a, n=160):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'

v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
print('gate atoms behind the two residuals:')
for a in (21616, 29538):
    print(f'  a{a}: {pr(a)}   defines x_{atom_out.get(a)}')
    for w in sorted(set(L.avars[a])):
        print(f'     x_{w:<7} {"FREE" if w in FREE else "gate":<5} '
              f'val={str(v0[w])[:26]:<28} bits={v0[w].bit_length():<5} '
              f'== p? {v0[w] == P}  consumers {len(L.var_atoms[w])}')

for target, gate, handles in [(21617, 21616, (986, 5040)), (29539, 29538, (11360, 30163))]:
    print(f'\n===== closing a{target} =====')
    av0 = L.all_atom_values(v0)
    for u, delta in newton_moves(target, v0, av0, topn=30)[:12]:
        v = list(v0); v[u] = v[u] + delta
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        if av[target] % P != 0: continue
        # required value of the gate's defined variable
        w = atom_out[gate][1]
        tgt = T.solve_lin(target, w, v)
        if tgt is None:
            print(f'  x_{u}: cannot solve a{target} for x_{w}'); continue
        vv = list(v); vv[w] = tgt
        done = False
        for h in handles:
            if h not in FREE: continue
            nv = T.solve_lin(gate, h, vv)
            if nv is None:
                print(f'  x_{u}: handle x_{h} does not divide '
                      f'(need x_{w}={str(tgt)[:20]}...)')
                continue
            tr = list(v); tr[h] = nv
            ad.fwd(tr, rounds=6)
            at = L.all_atom_values(tr)
            s = score(tr)
            nz = [b for b in range(L.NA) if at[b]]
            print(f'  x_{u} + handle x_{h}: score {s} ({s-score(v0):+d})  nonzero {nz}')
            if s > score(v0):
                T.save(tr, os.path.join(HERE, f'absorb_{s}.json'))
                print(f'     *** saved absorb_{s}.json')
            done = True
        if not done:
            print(f'  x_{u}: no usable handle among {handles}')
