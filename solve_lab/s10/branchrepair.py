"""S10 step 92: branch (1,1) + repair of the free inputs that branch pins."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    return s if c == 1 else ('-' + s if c == -1 else f'{c}*{s}')
def pretty(a):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')
    return o if len(o) < 190 else o[:190] + ' ...'

for a in [35758, 3568, 3570, 19088, 22233, 22235, 36602, 37887, 7930, 29539, 40826, 41512]:
    print(f'a{a:<6} eqs={len(L.atom2eq[a]):>3} '
          f'{"->x_"+str(atom_out[a]) if a in atom_out else "CHECK":>12}  {pretty(a)}')

base = L.load(os.path.join(HERE, 'forward_state.json'))
v = list(base); v[2081] = 1; v[4287] = 1
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
print(f'\nbranch(1,1) raw: failing {len(L.failing_eqs(av))}')

print('\n=== greedy free-input repair ===')
for it in range(14):
    av = L.all_atom_values(v)
    fail = L.failing_eqs(av)
    nz = [a for a in range(L.NA) if av[a]]
    moved = False
    for a in nz:
        for w in sorted(set(L.avars[a])):
            if w not in FREE: continue
            nv = T.solve_lin(a, w, v)
            if nv is None or nv == v[w]: continue
            trial = list(v); trial[w] = nv
            ad.fwd(trial, rounds=6)
            at = L.all_atom_values(trial)
            ft = L.failing_eqs(at)
            if len(ft) < len(fail):
                print(f'  it{it}: a{a} solved via FREE x_{w} -> failing '
                      f'{len(fail)} -> {len(ft)}  (score {L.NEQ-len(ft)})')
                v = trial; moved = True; break
        if moved: break
    if not moved:
        print(f'  it{it}: no single free-input move improves; '
              f'failing {len(fail)}, nonzero {nz}')
        break
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
print(f'\nfinal: failing {len(fail)}  score {L.NEQ-len(fail)}  '
      f'nonzero {[a for a in range(L.NA) if av[a]]}')
T.save(v, os.path.join(HERE, 'branch11.json'))
