"""Is there any direction that both moves the post-solve class AND preserves solvability?

Answered structurally, not by more sampling -- my own starvation table says sampling cannot
answer it.  Every lat3.analyse line already records the shape of the other-rows system:

    knobs=K  other-rows=M  ... FEASIBLE/INFEASIBLE  kernel-dim=d

from which rank(D_other) = K - d and the ROW DEFICIENCY is  M - (K - d).

Hypothesis: the other-rows system at cfg0 has FULL ROW RANK (deficiency 0), hence is solvable for
ANY right-hand side; leaving the affine span breaks atoms, which adds rows faster than it adds
knobs, driving the deficiency positive and making the system over-determined and generically
infeasible.  If that holds across every configuration ever logged, the infeasibility is intrinsic
to leaving the span rather than an artifact of the particular selectors relaxed.

Reads only my own run logs.
"""
import sys, re, glob, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')

pat = re.compile(r'knobs=(\d+) other-rows=(\d+) target-rows=\[[^\]]*\]: (FEASIBLE|INFEASIBLE) kernel-dim=(\d+)')
rows = []
for fn in sorted(glob.glob('runs_*.log')):
    for line in open(fn, errors='ignore'):
        m = pat.search(line)
        if m:
            K, M, feas, d = int(m.group(1)), int(m.group(2)), m.group(3), int(m.group(4))
            rank = K - d
            rows.append((fn, K, M, d, rank, M - rank, feas == 'FEASIBLE'))

print("configurations recovered from my own logs: %d\n" % len(rows), flush=True)
tab = collections.Counter()
for fn, K, M, d, rank, defi, feas in rows:
    tab[(defi, feas)] += 1

print("%-12s %-10s %s" % ("deficiency", "feasible", "count"), flush=True)
for (defi, feas), n in sorted(tab.items()):
    print("%-12d %-10s %d" % (defi, feas, n), flush=True)

zero_feas = sum(n for (defi, f), n in tab.items() if defi == 0 and f)
zero_inf = sum(n for (defi, f), n in tab.items() if defi == 0 and not f)
pos_feas = sum(n for (defi, f), n in tab.items() if defi > 0 and f)
pos_inf = sum(n for (defi, f), n in tab.items() if defi > 0 and not f)
print("\ndeficiency == 0 : feasible %d, infeasible %d" % (zero_feas, zero_inf), flush=True)
print("deficiency >  0 : feasible %d, infeasible %d" % (pos_feas, pos_inf), flush=True)
print("\nHYPOTHESIS 'deficiency == 0  <=>  feasible':", flush=True)
if zero_inf == 0 and pos_feas == 0:
    print("  HOLDS on every configuration logged (%d)." % len(rows), flush=True)
else:
    print("  FAILS: %d zero-deficiency infeasible, %d positive-deficiency feasible."
          % (zero_inf, pos_feas), flush=True)

# distinct system shapes
shapes = collections.Counter((K, M, d) for _, K, M, d, _, _, _ in rows)
print("\ndistinct (knobs, other-rows, kernel-dim) shapes seen: %d" % len(shapes), flush=True)
for (K, M, d), n in shapes.most_common(12):
    print("   knobs=%-4d rows=%-4d kerdim=%-3d rank=%-4d deficiency=%-3d  x%d"
          % (K, M, d, K - d, M - (K - d), n), flush=True)
