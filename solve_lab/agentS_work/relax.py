"""Relaxed selectors: the remaining way to generate VALID test cases.

Why this and not the trade knobs: trade.py confirms all four trade knobs lie inside the 54-knob
affine span lat3.analyse already optimises over, so displacing along them cannot change the
membership answer.  Selectors are different -- they are non-affine (saturating), so taking one
OFF {0,1} changes the measured structure itself, which is the only thing that can move the answer.

Lead from agent R (unadjudicated, another model's parse): relaxing a selector off {0,1} does not
force its mux atoms nonzero -- only the boolean-ness atoms are forced.  Tested here in MY parse
rather than taken on faith; R's atom indices are not comparable to mine so nothing is imported.

PRIMARY REPORT: count of VALID cases = other rows re-solvable AND the target mod-p class moved.
"""
import sys, json, collections, pickle, time, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat2, lat3
import harness as H, engine as E, fast
P = C.P
TGT = [20215, 28647]

seed0 = dict(C.BASE)
res = lat3.analyse(seed0, 'cfg0')
y, ker, n0, knobs, aff, v0, bad0, atoms = res
base = dict(seed0)
for j, f in enumerate(knobs):
    if n0[j]:
        base[f] = v0[f] + n0[j]
vb = E.forward(base); badb = E.badatoms(vb)
nsb = {'v': vb, '__builtins__': {}}
class0 = tuple(eval(H.acodes[a], nsb) % P for a in TGT)
print("base near-solution: bad=%s fails=%d" % (sorted(badb), len(E.eqfails(badb))), flush=True)

SEL = sorted(f for f in C.cluster_cone() if C.isbool(f))
random.seed(31337)
sample = random.sample(SEL, 14)
vals = [2, -1, 5]

valid = blocked = solved = 0
infeas = sameclass = 0
best = (10**9, None)
print("\n=== relaxing %d selectors x %d off-boolean values ===" % (len(sample), len(vals)), flush=True)
for s in sample:
    for val in vals:
        ns = dict(base); ns[s] = val
        try:
            v = E.forward(ns); bad = E.badatoms(v); ff = E.eqfails(bad)
            nsx = {'v': v, '__builtins__': {}}
            cls = tuple(eval(H.acodes[a], nsx) % P for a in TGT)
        except Exception:
            continue
        moved = (cls != class0)
        newbad = sorted(set(bad) - set(badb))
        if len(ff) < best[0]:
            best = (len(ff), (s, val))
        try:
            r = lat3.analyse(ns, '  x_%d=%d' % (s, val))
        except Exception as e:
            print("  x_%d=%-3d analyse ERR %s" % (s, val, type(e).__name__), flush=True); continue
        if r is None:
            infeas += 1
            print("  x_%d=%-3d moved=%-5s newbad=%s fails=%d -> other rows INFEASIBLE"
                  % (s, val, moved, newbad[:6], len(ff)), flush=True)
            continue
        if not moved:
            sameclass += 1
            print("  x_%d=%-3d solvable but CLASS UNMOVED (not valid)" % (s, val), flush=True)
            continue
        valid += 1
        if r[0] is not None:
            solved += 1
            print("  x_%d=%-3d *** VALID and SOLVED -- obstruction cfg0-local ***" % (s, val), flush=True)
            json.dump({str(x): str(int(z)) for x, z in ns.items()}, open('S_relax_hit.json', 'w'))
        else:
            blocked += 1
            print("  x_%d=%-3d VALID, blocked (fails=%d)" % (s, val, len(ff)), flush=True)

print("\n=== SUMMARY (relaxed selectors) ===", flush=True)
print("attempts                   : %d" % (len(sample) * len(vals)), flush=True)
print("other rows infeasible      : %d" % infeas, flush=True)
print("solvable but class unmoved : %d" % sameclass, flush=True)
print("VALID CASES                : %d   (blocked=%d solved=%d)" % (valid, blocked, solved), flush=True)
print("best fails seen            : %s at %s" % (best[0], best[1]), flush=True)
if valid == 0:
    print("=> relaxed selectors starve too.", flush=True)
elif solved == 0:
    print("=> %d valid cases, all blocked." % valid, flush=True)
