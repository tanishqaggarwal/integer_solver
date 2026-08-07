"""P's test: how many atoms touch two or more distinct selector variables?

Measured TWO ways, because that is where I expect the decomposition difference to live:
  (a) ATOM-LOCAL   -- selector variables appearing literally in the atom's own expression.
                      This is what "atoms containing two or more selector variables" reads as.
  (b) COMPOSED     -- selectors appearing in the atom's transitive CONE (i.e. reaching the atom
                      through chains of defined intermediate variables).
My saturation measurement ranges over composed residual values, not over atom-local text, so if
(a) is 0 and (b) is large the two models are describing the same object at different depths.
"""
import sys, re, collections, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H, engine as E

SEL = set(f for f in C.cluster_cone() if C.isbool(f))
print("selector set (my parse): %d boolean selectors in the cluster cone" % len(SEL), flush=True)

# ---------- (a) atom-local ----------
local = collections.Counter()
examples = []
occ = collections.Counter()
for a, t in enumerate(H.atoms):
    if not t:
        continue
    vs = {int(m.group(1)) for m in re.finditer(r'x_(\d+)', t)}
    s = vs & SEL
    for v in s:
        occ[v] += 1
    if len(s) >= 2:
        local[len(s)] += 1
        if len(examples) < 8:
            examples.append((a, sorted(s), t[:140]))
print("\n(a) ATOM-LOCAL: atoms containing >=2 distinct selectors: %d" % sum(local.values()), flush=True)
if examples:
    for a, s, t in examples:
        print("     a%-6d selectors=%s  %s" % (a, s, t), flush=True)
print("    atoms-per-selector (atom-local occurrence count): min=%d max=%d mean=%.2f"
      % (min(occ.values()), max(occ.values()), sum(occ.values()) / len(occ)), flush=True)
print("    distribution:", dict(sorted(collections.Counter(occ.values()).items())), flush=True)

# ---------- (b) composed ----------
print("\n(b) COMPOSED: selectors reaching each cluster row through its transitive cone", flush=True)
for r in C.ROWS:
    order, fr, seen = E.cone(r)
    s = set(fr) & SEL
    print("    a%-6d cone free vars=%-5d of which selectors=%d" % (r, len(fr), len(s)), flush=True)

# how many atoms have >=2 selectors in their cone?
multi = 0
tot = 0
hist = collections.Counter()
for a in range(len(H.atoms)):
    if not H.atoms[a]:
        continue
    try:
        order, fr, seen = E.cone(a)
    except Exception:
        continue
    tot += 1
    n = len(set(fr) & SEL)
    hist[min(n, 10)] += 1
    if n >= 2:
        multi += 1
    if tot >= 4000:
        break
print("    of %d atoms sampled: %d have >=2 selectors in their cone" % (tot, multi), flush=True)
print("    histogram of #selectors-in-cone (10 = 10 or more):", dict(sorted(hist.items())), flush=True)
