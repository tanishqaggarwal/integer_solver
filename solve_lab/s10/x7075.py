import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    return s if c == 1 else ('-' + s if c == -1 else f'{c}*{s}')
def pretty(a):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')
    return o if len(o) < 200 else o[:200] + ' ...'

# how many variables carry each "wire" constant?
cnt = collections.Counter(v)
print('most common variable VALUES across all 38,748 variables:')
for val, n in cnt.most_common(8):
    print(f'   value {str(val)[:34]:<36} carried by {n} variables  (bits {val.bit_length()})')

print('\n=== the x_7075 chain ===')
seen, stack = set(), [7075]
while stack:
    t = stack.pop()
    if t in seen: continue
    seen.add(t)
    a = definer.get(t)
    print(f'  x_{t:<7} = {str(v[t])[:30]:<32} '
          f'{"def by a"+str(a) if a is not None else "FREE INPUT"}'
          f'{"  : " + pretty(a) if a is not None else ""}')
    if a is None: continue
    for w in L.avars[a]:
        if w != t and w not in seen: stack.append(w)
    if len(seen) > 25: break
print(f'\nx_7075 appears in atoms: {sorted(L.var_atoms[7075])}')
for a in sorted(L.var_atoms[7075]):
    print(f'   a{a:<6} eqs={len(L.atom2eq[a]):>3} '
          f'{"->x_"+str(atom_out[a]) if a in atom_out else "CHECK"}  {pretty(a)}')
