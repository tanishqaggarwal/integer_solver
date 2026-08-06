"""WR step 2: drive the whole wire to w inside frame3+root and measure."""
import os, sys, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = W.wire_of(base)
F = W.F_WIRE
print(f'wire {len(WIRE)}; detached vars in wire: '
      f'{[u for u in F.detach if u in set(WIRE)]}')

b2 = list(base); F.fwd(b2)
F.report(b2, 'baseline')

vals = []
for x in sys.argv[1:]:
    vals.append(eval(x, {'P': P}))
if not vals:
    vals = [1, -1, 2, -2, 3, 0, P, -P, 2 * P, P * P, P + 1, P - 1]

out = {}
for Wv in vals:
    v = list(b2)
    for u in WIRE:
        v[u] = Wv
    F.fwd(v, rounds=10)
    held = sum(1 for u in WIRE if v[u] == Wv)
    av, nz, fail, sc = F.report(v, f'w={str(Wv)[:30]} held={held}/{len(WIRE)}')
    out[str(Wv)] = sc
    tag = str(Wv)
    if len(tag) > 12:
        tag = 'big' + str(abs(hash(tag)) % 10000)
    T.save(v, os.path.join(HERE, 'wr_w_%s.json' % tag.replace('-', 'm')))
print()
for k, s in out.items():
    print(f'   w={k[:40]:<42} {s}')
