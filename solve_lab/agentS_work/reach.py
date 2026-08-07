"""Is my BFS closure LOCAL or GLOBAL?

bfs.py explored from cfg0 by single flips and closed at 48 mod-p 5-tuples.  That is exhaustive
over the closure reachable from cfg0 under single flips.  It is NOT obviously exhaustive over the
2^256 selector domain.  This script evaluates configurations the BFS could not have reached --
random subsets of every weight, each one FULLY specified (every one of the 256 selectors set
explicitly, not a few flips off cfg0) -- and asks whether the resulting tuple is one of the 48.

Uses bfs.py's own key_of so the tuples are directly comparable.
"""
import sys, json, collections, pickle, random, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H, engine as E, fast
P = C.P

base = dict(C.BASE)
v0 = E.forward(base)
ROWS = C.ROWS
BOOLS = sorted(f for f in C.cluster_cone() if C.isbool(f))
print("selectors: %d ; rows: %s" % (len(BOOLS), ROWS), flush=True)

seen = pickle.load(open('bfs_image.pkl', 'rb'))
KNOWN = set(seen.keys())
print("BFS image (closed by exhaustion under single flips): %d tuples" % len(KNOWN), flush=True)


def key_of(assign):
    v, _ = fast.apply_delta(v0, assign)
    ns = {'v': v, '__builtins__': {}}
    return tuple(eval(H.acodes[a], ns) % P for a in ROWS)


random.seed(20260807)
sizes = [1, 2, 5, 17, 32, 64, 128, 192, 200, 256]
TRIALS = 30
newtuples = {}
stats = []
a20215_vals = collections.Counter()
t0 = time.time()
for sz in sizes:
    inside = outside = 0
    for t in range(TRIALS):
        S = set(random.sample(BOOLS, sz))
        assign = {f: (1 if f in S else 0) for f in BOOLS}
        # also randomise the two switches
        assign[30163] = random.choice([0, 1])
        assign[11559] = random.choice([0, 1])
        try:
            k = key_of(assign)
        except Exception:
            continue
        a20215_vals[k[3]] += 1
        if k in KNOWN:
            inside += 1
        else:
            outside += 1
            if k not in newtuples:
                newtuples[k] = (sz, sorted(S)[:6], assign[30163], assign[11559])
    stats.append((sz, inside, outside))
    print("  |S|=%-4d  inside the 48: %-3d   OUTSIDE: %-3d   (new tuples so far %d)  %.0fs"
          % (sz, inside, outside, len(newtuples), time.time() - t0), flush=True)

print("\n=== VERDICT ===", flush=True)
tot_in = sum(s[1] for s in stats)
tot_out = sum(s[2] for s in stats)
print("total trials %d : inside %d, outside %d" % (tot_in + tot_out, tot_in, tot_out), flush=True)
if tot_out == 0:
    print("ALL sampled configurations landed inside the 48-tuple image.", flush=True)
    print("=> the closure is GLOBAL on this evidence; sec 3 binds the instance.", flush=True)
else:
    print("%d NEW tuples outside the BFS closure => the closure is LOCAL, not global." % len(newtuples), flush=True)
    print("=> sec 3's 'a20215 never 0' does NOT bind the instance.", flush=True)

print("\ndistinct a20215 values mod p seen across ALL random trials: %d" % len(a20215_vals), flush=True)
print("  a20215 == 0 ever?  %s" % (0 in a20215_vals), flush=True)
for v, n in a20215_vals.most_common(8):
    print("    %s   x%d" % (str(v)[:50], n), flush=True)
pickle.dump({'new': newtuples, 'stats': stats, 'a20215': dict(a20215_vals)}, open('reach.pkl', 'wb'))
