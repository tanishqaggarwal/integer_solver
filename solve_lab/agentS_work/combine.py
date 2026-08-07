"""Combine kernel2 (img0-13) and dirsearch (img13-47) into the count that matters.

INDEPENDENT TEST CASE = other rows solvable AND post-solve target class distinct.
The post-solve class is the criterion; the pre-solve class is the trap and is never used.
"""
import sys, re, glob, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')

cases = {}   # img index -> (status, class tuple)

# ---- dirsearch RESULT lines (machine readable) ----
for fn in glob.glob('runs_dir_*.log'):
    for line in open(fn, errors='ignore'):
        m = re.match(r'RESULT (\d+) (\S+) (\S+) (\S+)', line)
        if not m:
            continue
        i, status, c1, c2 = int(m.group(1)), m.group(2), m.group(3), m.group(4)
        if status == 'other-rows-infeasible':
            cases[i] = ('infeasible', None)
        elif status in ('blocked', 'SOLVED'):
            cases[i] = (status, (c1, c2))

# ---- kernel2 log (block aware) ----
lines = open('runs_kernel2.log', errors='ignore').read().split('\n')
for idx, line in enumerate(lines):
    m = re.search(r'\[\s*img(\d+)\(\|on\|=\d+\)\] knobs=\d+ other-rows=\d+ target-rows=\[[^\]]*\]: (FEASIBLE|INFEASIBLE)', line)
    if not m:
        continue
    i, feas = int(m.group(1)), m.group(2)
    if i in cases:
        continue
    if feas == 'INFEASIBLE':
        cases[i] = ('infeasible', None); continue
    c1 = c2 = None
    solved = False
    for j in range(idx, min(idx + 26, len(lines))):
        if 'YES -> FULL SOLVE' in lines[j]:
            solved = True
        a = re.search(r'row a20215: residual mod p = (\d+)', lines[j])
        b = re.search(r'row a28647: residual mod p = (\d+)', lines[j])
        if a: c1 = a.group(1)
        if b: c2 = b.group(1)
        if re.search(r'\[\s*img\d+', lines[j]) and j > idx:
            break
    cases[i] = ('SOLVED' if solved else 'blocked', (c1, c2) if c1 else None)

tot = len(cases)
infeas = [i for i, (s, c) in cases.items() if s == 'infeasible']
feas = [i for i, (s, c) in cases.items() if s in ('blocked', 'SOLVED')]
solved = [i for i, (s, c) in cases.items() if s == 'SOLVED']
blocked = [i for i, (s, c) in cases.items() if s == 'blocked']

classes = collections.OrderedDict()
for i in sorted(feas):
    c = cases[i][1]
    if c is None:
        continue
    classes.setdefault(c, []).append(i)

print("=== DEFICIENCY-DIRECTED SWEEP: combined result ===\n", flush=True)
print("image points analysed          : %d of 48" % tot, flush=True)
print("other rows infeasible          : %d  (not test cases at all)" % len(infeas), flush=True)
print("other rows solvable            : %d" % len(feas), flush=True)
print("   of those, blocked           : %d" % len(blocked), flush=True)
print("   of those, SOLVED            : %d" % len(solved), flush=True)
print("\nINDEPENDENT test cases (distinct POST-SOLVE class): %d" % len(classes), flush=True)
for c, idxs in classes.items():
    st = set(cases[i][0] for i in idxs)
    print("   class a20215=%s... a28647=%s...  imgs=%s  %s"
          % (str(c[0])[:26], str(c[1])[:26], idxs, '/'.join(sorted(st))), flush=True)
print("\nstarvation rate (infeasible / analysed): %.0f%%" % (100.0 * len(infeas) / max(tot, 1)), flush=True)
if solved:
    print("\n=> SOLVED at %s -- the endgame condition DISSOLVES." % solved, flush=True)
elif len(classes) >= 5:
    print("\n=> %d independent classes, ALL BLOCKED: configuration-independence of the joint p*Z^2"
          % len(classes), flush=True)
    print("   obstruction, at the resolution this sweep reaches.", flush=True)
else:
    print("\n=> only %d independent classes -- directed search STARVED too; not enough to claim"
          % len(classes), flush=True)
    print("   configuration-independence. Report the rate, not a conclusion.", flush=True)
