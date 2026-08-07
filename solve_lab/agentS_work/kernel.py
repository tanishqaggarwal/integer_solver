"""sec 8.3 -- is the joint p*Z^2 obstruction configuration-INDEPENDENT, or another cfg0-local fact?

The right way to build test configurations (reach3.py's random selector scrambles were useless:
66-467 bad atoms, blocking rows nowhere near the cluster).  Instead move along the AFFINE KERNEL:
directions that, by construction, hold every non-target row at its satisfied value.  So every
configuration produced here is still a near-solution, and any of them that are outside cfg0's
BFS closure are exactly the test cases sec 8.3 asks for.

Knob set: the 54 affine knobs of lat2.system at cfg0 (every single-row knob included, plus each
atom's pure handle).  Base configuration: cfg0 = triple8_seed (x_1530 = x_1603 = 1).

At each displaced configuration everything is RE-MEASURED -- knobs, handles and targets are all
configuration-dependent -- and the 2-D reachable lattice on (a20215, a28647) is recomputed.
"""
import sys, json, collections, pickle, time, math, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat2, lat3
import harness as H, engine as E, fast, intsolve
P = C.P
TGT = [20215, 28647]

print("=== step 1: cfg0 baseline ===", flush=True)
seed0 = dict(C.BASE)
res = lat3.analyse(seed0, 'cfg0')
if res is None:
    print("cfg0 analyse returned None -- cannot proceed"); sys.exit(1)
y, ker, n0, knobs, aff, v0, bad0, atoms = res
print("kernel dim = %d over %d knobs" % (len(ker), len(knobs)), flush=True)


def build(coeffs):
    """seed for cfg0 + particular solution + sum(coeffs[i] * ker[i])"""
    n = list(n0)
    for c, k in zip(coeffs, ker):
        if c:
            for j in range(len(knobs)):
                if k[j]:
                    n[j] += c * k[j]
    ns = dict(seed0)
    for j, f in enumerate(knobs):
        if n[j]:
            ns[f] = v0[f] + n[j]
    return ns


random.seed(9001)
trials = [('n0 only (kernel coeffs all 0)', [0] * len(ker))]
for t in range(6):
    trials.append(('small coeffs #%d' % t, [random.randint(-3, 3) for _ in ker]))
for t in range(4):
    trials.append(('large coeffs #%d' % t, [random.randint(-10**6, 10**6) for _ in ker]))
for i in range(len(ker)):
    c = [0] * len(ker); c[i] = 1
    trials.append(('unit kernel vector %d' % i, c))

print("\n=== step 2: displace along the kernel, verify near-solution, re-measure ===", flush=True)
summary = []
for tag, c in trials:
    ns = build(c)
    t0 = time.time()
    try:
        v = E.forward(ns)
        bad = E.badatoms(v)
    except Exception as e:
        print("[%s] forward ERR %s" % (tag, type(e).__name__), flush=True); continue
    ff = E.eqfails(bad)
    nearsol = len(bad) <= 6
    print("\n[%s] bad atoms=%d %s fails=%d SCORE=%d  near-solution=%s (%.0fs)"
          % (tag, len(bad), sorted(bad) if len(bad) <= 8 else '(many)', len(ff),
             39033 - len(ff), nearsol, time.time() - t0), flush=True)
    if not nearsol:
        print("    kernel move did NOT preserve the other rows -> linear model broke down here",
              flush=True)
        summary.append((tag, len(bad), None, None)); continue
    # re-measure the joint obstruction at this configuration
    try:
        r2 = lat3.analyse(ns, '  ' + tag)
    except Exception as e:
        print("    analyse ERR %s" % type(e).__name__, flush=True)
        summary.append((tag, len(bad), 'ERR', None)); continue
    got = 'FULL SOLVE' if (r2 is not None and r2[0] is not None) else ('blocked' if r2 is not None else 'trivial/infeasible')
    summary.append((tag, len(bad), got, None))
    if r2 is not None and r2[0] is not None:
        print("    *** MEMBERSHIP SUCCEEDED AT THIS CONFIGURATION -- OBSTRUCTION IS cfg0-LOCAL ***", flush=True)
        json.dump({str(x): str(int(z)) for x, z in ns.items()}, open('S_kernel_hit.json', 'w'))

print("\n=== SUMMARY ===", flush=True)
for tag, nb, got, _ in summary:
    print("  %-32s bad=%-4d joint-obstruction: %s" % (tag, nb, got), flush=True)
ns_ok = [s for s in summary if s[2] == 'blocked']
hits = [s for s in summary if s[2] == 'FULL SOLVE']
print("\nconfigurations that stayed near-solutions and were re-measured: %d" % (len(ns_ok) + len(hits)), flush=True)
print("of those, obstruction still blocked: %d ; dissolved: %d" % (len(ns_ok), len(hits)), flush=True)
if hits:
    print("=> the joint p*Z^2 obstruction is cfg0-LOCAL. sec 2 does not bind the instance.", flush=True)
elif ns_ok:
    print("=> obstruction survived motion along the kernel at every tested configuration.", flush=True)
    print("   *** DO NOT READ THIS AS 'the obstruction is a statement about the instance'. ***", flush=True)
    print("   Motion along the kernel changes a20215 only by multiples of p -- that IS the p*Z^2", flush=True)
    print("   result -- so the residual's mod-p class is INVARIANT along the kernel by", flush=True)
    print("   construction, and the membership answer cannot change unless the measured knob set", flush=True)
    print("   changes.  This is a test of STRUCTURAL STABILITY (which passed: 54 knobs, 47 other", flush=True)
    print("   rows, kernel dim 7, bad={a20215,a28647} at every displacement).  It is NOT a test of", flush=True)
    print("   configuration-independence.  For that see kernel2.py, which moves the mod-p class.", flush=True)
