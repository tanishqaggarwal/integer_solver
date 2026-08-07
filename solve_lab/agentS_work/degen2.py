"""Degeneracy discriminator.

P's degenerate family = a merge sees two equal live inputs, so A = B = 0 and BOTH of that
block's congruences vanish IDENTICALLY.  The observable signature is not "the row is zero"
(that is just the congruence being satisfied) but "the row is zero AND has stopped responding
to every knob in its cone" -- an identically-vanishing constraint cannot be moved.

So: over my converged BFS image, for every configuration and every cluster row, measure
  (value == 0)  and  (number of knobs that still move the row).
A degenerate block shows value 0 with responsiveness 0.
"""
import sys, json, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H, engine as E, fast
P = C.P

seen = pickle.load(open('bfs_image.pkl', 'rb'))
cone = C.cluster_cone()
probe = sorted(cone)
print("BFS image configurations: %d ; cone knobs probed per config: %d" % (len(seen), len(probe)), flush=True)

rows = C.ROWS
tally = collections.Counter()
frozen_cases = []
for i, (k, a) in enumerate(sorted(seen.items(), key=lambda kv: len(kv[1]))):
    seed = dict(C.BASE); seed.update(a)
    v0 = E.forward(seed)
    ns = {'v': v0, '__builtins__': {}}
    val = {r: eval(H.acodes[r], ns) for r in rows}
    # responsiveness: how many cone knobs change this row at all
    resp = collections.Counter()
    for f in probe:
        try:
            v1, _ = fast.apply_delta(v0, {f: (0 if v0[f] else 1) if C.isbool(f) else v0[f] + 1})
        except Exception:
            continue
        ns1 = {'v': v1, '__builtins__': {}}
        for r in rows:
            try:
                if eval(H.acodes[r], ns1) != val[r]:
                    resp[r] += 1
            except Exception:
                pass
    for r in rows:
        z = (val[r] == 0)
        fr = (resp[r] == 0)
        tally[(r, z, fr)] += 1
        if z and fr:
            frozen_cases.append((i, r, len(a)))
    if i < 6 or (i % 12 == 0):
        print("  cfg%-3d zero=%s  responsiveness=%s"
              % (i, [r for r in rows if val[r] == 0], {("a%d" % r): resp[r] for r in rows}), flush=True)

print("\n=== summary over %d configurations ===" % len(seen), flush=True)
for r in rows:
    print("  a%-6d exactly0&frozen(DEGENERATE)=%d  exactly0&responsive=%d  nonzero&responsive=%d  nonzero&frozen=%d"
          % (r, tally[(r, True, True)], tally[(r, True, False)],
             tally[(r, False, False)], tally[(r, False, True)]), flush=True)
print("\nDEGENERATE (zero and frozen) cases found:", len(frozen_cases), flush=True)
if frozen_cases:
    print("  ", frozen_cases[:20], flush=True)
