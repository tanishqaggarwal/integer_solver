import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from repair2 import candidates, score, FORBID
P = ad.P
definer, atom_out = L.definer, L.atom_out
def term(m, c):
    if not m: return f'{c}'
    s = '*'.join(f'x_{w}' for w in m)
    return s if c == 1 else ('-' + s if c == -1 else f'{c}*{s}')
def pretty(a, n=170):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(term(m, c) for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'

v = L.load(os.path.join(HERE, 'br10.json'))
s0 = score(v)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'br10 score {s0}; nonzero {nz}')
for a in nz:
    print(f'\na{a} ({len(L.atom2eq[a])} eqs) {pretty(a)}')
    cs = candidates(a, v)
    if not cs: print('   NO candidate move'); continue
    for mv in cs:
        trial = list(v)
        for u, uv in mv: trial[u] = uv
        ad.fwd(trial, rounds=6)
        s = score(trial)
        at = L.all_atom_values(trial)
        lbl = ','.join(f'x_{u}' for u, _ in mv)
        nzz = [b for b in range(L.NA) if at[b]]
        print(f'   move {lbl:<18} -> score {s} ({s-s0:+d})  nonzero {nzz}')
