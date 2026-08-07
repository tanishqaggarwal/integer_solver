"""W stage 17: block 7181's six slots are ALL free inputs and its gate is DEAD.
Its outputs are nonzero, which breaks its two off-pins.  Does zeroing them help?"""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st0 = frameB.State(fr, fv)
print('base score', st0.score(), 'fails', sorted(st0.fails))
CAND = [9118, 8731, 29854, 31864, 642, 28730, 6418, 12553, 31861, 14865]
print()
for u in CAND:
    print('  x_%-6d free=%s  val=%s' % (u, u in fr.free, str(v0[u])[:40]))
print()
best = (st0.score(), {})
for r in range(1, 5):
    for comb in itertools.combinations(CAND, r):
        s = st0.clone().set_free({u: 0 for u in comb})
        if s.score() > best[0]:
            best = (s.score(), {u: 0 for u in comb})
            print('  IMPROVE', s.score(), comb, sorted(s.fails))
        elif s.score() >= st0.score():
            print('  equal  ', s.score(), comb, sorted(s.fails))
    print(' r=%d done, best=%d' % (r, best[0]), flush=True)
print('BEST', best[0])
