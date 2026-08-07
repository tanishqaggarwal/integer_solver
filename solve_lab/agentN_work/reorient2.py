"""Re-orientation, executed rather than argued.

Orienting a check atom `x_v - rest` into a DEFINITION of x_v forces that atom to 0 for every choice
of free inputs.  But inside a frame where x_v is already free, "force the atom to 0" and "choose the
value of x_v that makes the atom 0" are the same set of assignments -- and the second needs no new
frame.  So each legal re-orientation of a region atom is testable directly as a knob setting.

This finds, for every nonzero atom at the witness, a free knob with unit response, applies the
setting that zeroes that atom, and reports the ACTUAL score.  Then it tries the combinations.
"""
import json, itertools, os
import ev, model
import optN
from optN import make, POOL, WIT, FREE, fr, inner, atom_eqs
from orient import unit_targets

d = model.get()
atom_src = d['atom_src']
NZATOMS = [22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887]

st = make(WIT)
print('witness: score %d, failing %s' % (st.score(), sorted(st.fails)))
print('\natom values at the witness and the free knobs with unit response:')
moves = {}
for a in NZATOMS:
    val = st.av.get(a, 0)
    u = unit_targets(atom_src[a])
    entries = []
    for v, s in u.items():
        if v not in FREE:
            entries.append('x_%d NOT a frame free input' % v)
            continue
        h = st.clone().set_free({v: st.fv.get(v, 0) + 1})
        dd = h.av.get(a, 0) - val
        entries.append('x_%d (in frame, d(atom)/d = %+d)' % (v, dd))
        if dd in (1, -1):
            moves.setdefault(a, []).append((v, dd))
    print('  atom %-6d value %-14s targets: %s'
          % (a, ('0' if val == 0 else '%d digits' % len(str(abs(val)))),
             '; '.join(entries) if entries else 'NONE'))

print('\n=== apply each zeroing move on its own ===')
res = []
for a, lst in moves.items():
    for v, dd in lst:
        val = st.av.get(a, 0)
        if val == 0:
            continue
        newv = st.fv.get(v, 0) - val // dd
        if (val % dd) != 0:
            print('  atom %d via x_%d: not an integer step' % (a, v))
            continue
        g = st.clone().set_free({v: newv})
        ok = g.av.get(a, 0) == 0
        print('  zero atom %-6d via x_%-6d -> atom now %s, score %d (was %d), failing %d'
              % (a, v, 'ZERO' if ok else 'NONZERO', g.score(), st.score(), len(g.fails)))
        res.append(dict(atom=a, knob=v, zeroed=ok, score=g.score()))

print('\n=== all combinations of the independent zeroing moves ===')
flat = []
for a, lst in moves.items():
    if st.av.get(a, 0) != 0:
        flat.append((a, lst[0][0], lst[0][1]))
print('available moves: %s' % [(a, 'x_%d' % v) for a, v, _ in flat])
best = (st.score(), ())
for r in range(1, len(flat) + 1):
    for C in itertools.combinations(flat, r):
        g = st.clone()
        okall = True
        for a, v, dd in C:
            val = g.av.get(a, 0)
            if val % dd:
                okall = False
                break
            g = g.set_free({v: g.fv.get(v, 0) - val // dd})
        if not okall:
            continue
        s = g.score()
        zs = [a for a, v, dd in C if g.av.get(a, 0) == 0]
        if s > best[0]:
            best = (s, tuple((a, v) for a, v, _ in C))
        print('  %-46s -> score %d  (atoms actually zeroed: %s)'
              % (str([(a, 'x_%d' % v) for a, v, _ in C]), s, zs))
print('\nbest over all re-orientation combinations: score %d via %s' % best)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 'runs', 'reorient2.json'), 'w'), indent=1)
