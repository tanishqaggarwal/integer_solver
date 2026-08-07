"""Which logged configurations LEFT the span (shape != cfg0's 54/47/7) and still had full row
   rank (deficiency 0)?  Those are precisely the candidate directions that could both move the
   post-solve class and preserve solvability -- the bounded structural question.

   Reads only my own run logs; then re-derives the post-solve class for any it finds.
"""
import sys, re, glob, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')

pat = re.compile(r'\[([^\]]*)\] knobs=(\d+) other-rows=(\d+) target-rows=\[[^\]]*\]: (FEASIBLE|INFEASIBLE) kernel-dim=(\d+)')
CFG0 = (54, 47, 7)
found = []
for fn in sorted(glob.glob('runs_*.log')):
    lines = open(fn, errors='ignore').read().split('\n')
    for i, line in enumerate(lines):
        m = pat.search(line)
        if not m:
            continue
        label, K, M, feas, d = m.group(1).strip(), int(m.group(2)), int(m.group(3)), m.group(4), int(m.group(5))
        rank = K - d
        defi = M - rank
        if (K, M, d) == CFG0:
            continue
        # post-solve class, if the analyse block printed it
        cls = None
        for j in range(i, min(i + 24, len(lines))):
            mm = re.search(r'row a20215: residual mod p = (\d+)', lines[j])
            if mm:
                cls = mm.group(1); break
        found.append((fn, label, K, M, d, rank, defi, feas, cls))

print("logged configurations that LEFT cfg0's shape: %d\n" % len(found), flush=True)
print("%-22s %-26s %5s %5s %4s %5s %5s %-11s %s"
      % ("log", "label", "knobs", "rows", "ker", "rank", "defic", "feasible", "post-solve a20215 class"), flush=True)
zero_def = []
for fn, label, K, M, d, rank, defi, feas, cls in found:
    mark = ' <== FULL ROW RANK' if defi == 0 else ''
    print("%-22s %-26s %5d %5d %4d %5d %5d %-11s %s%s"
          % (fn[:22], label[:26], K, M, d, rank, defi, feas,
             (cls[:24] + '...') if cls else '-', mark), flush=True)
    if defi == 0:
        zero_def.append((fn, label, K, M, d, feas, cls))

print("\n=== ANSWER ===", flush=True)
print("configurations outside cfg0's shape WITH full row rank (deficiency 0): %d" % len(zero_def), flush=True)
CFG0CLASS = '22981624690591324143788809642515852940280603493270692712106986169263210356252'
for fn, label, K, M, d, feas, cls in zero_def:
    same = (cls == CFG0CLASS) if cls else None
    print("   [%s] %s  shape=(%d,%d,%d) feasible=%s  class==cfg0's? %s"
          % (fn, label, K, M, d, feas, same), flush=True)
if not zero_def:
    print("   NONE. Every direction that left the span lost full row rank, and every", flush=True)
    print("   positive-deficiency system in the logs was infeasible.", flush=True)
