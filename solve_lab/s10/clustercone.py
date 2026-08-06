"""S10 step 113: the ancestor cone of the CLUSTER's computed sides.

a21617 needs x_27522 == x_14623 (mod p);  a29539 needs x_1308 == x_14853 (mod p).
Reading the 29-variable cone is what cracked D0/K2 open.  Do the same here.
"""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
def pr(a, n=120):
    ts = sorted(L.polys[a].items(), key=lambda kv: (len(kv[0]), kv[0]))
    o = ' + '.join(('*'.join(f'x_{z}' for z in m) if c == 1 else
                    ('-' + '*'.join(f'x_{z}' for z in m) if c == -1 else
                     f'{c}*' + '*'.join(f'x_{z}' for z in m)) if m else str(c))
                   for m, c in ts).replace('+ -', '- ')
    return o if len(o) < n else o[:n] + ' ...'

for tgt, name in [(27522, 'a21617'), (1308, 'a29539'), (6858, 'a33796'),
                  (25442, 'a7930')]:
    cone, stack = set(), [tgt]
    while stack:
        t = stack.pop()
        if t in cone: continue
        cone.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: stack.append(w)
    fr = sorted(w for w in cone if w in FREE)
    print(f'\n=== cone of x_{tgt} ({name}): {len(cone)} vars, {len(fr)} free ===')
    if len(cone) <= 40:
        for t in sorted(cone):
            a = definer.get(t)
            tag = 'FREE' if a is None else f'a{a}'
            print(f'   x_{t:<7} {tag:<8} {str(v[t])[:24]:<26} bits={v[t].bit_length():<5}'
                  f'{"  : " + pr(a) if a is not None else ""}')
    else:
        print(f'   free inputs: {fr[:24]}{" ..." if len(fr) > 24 else ""}')
        wire = [t for t in cone if v[t] == P]
        print(f'   wire members (== p) in cone: {len(wire)}')
        print(f'   depth-1: {sorted(set(L.avars[definer[tgt]]) - {tgt})}')
