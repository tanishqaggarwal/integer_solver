"""Decisive test for the tension with Q's existence result.

p is PRIME.  An affine knob's effect on a row scales linearly with its step, so if ANY affine
knob f moves a20215 by an amount d with d != 0 (mod p), then n*d ranges over ALL of Z_p as n
ranges over Z -- and a20215 == 0 (mod p) is reachable, which would resolve the tension in Q's
favour and refute my sec 3 outright.

Conversely if every knob in the whole cone moves a20215 by a multiple of p, then a20215 mod p is
an INVARIANT of the affine directions and only the (saturating) selectors can change it -- which
is what my sec 2 lattice computation found at cfg0, here re-tested at configurations that my BFS
demonstrably could not reach.
"""
import sys, json, collections, pickle, random, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C
import harness as H, engine as E, fast
P = C.P

BOOLS = sorted(f for f in C.cluster_cone() if C.isbool(f))
CONE = sorted(C.cluster_cone())
SWITCH = {30163, 11559}
AFFINE = [f for f in CONE if f not in SWITCH and not C.isbool(f)]
print("cone knobs %d = %d selectors + %d switches + %d affine" %
      (len(CONE), len(BOOLS), len(SWITCH & set(CONE)), len(AFFINE)), flush=True)
TARGET = 20215

random.seed(11)
configs = [('cfg0', {})]
for i, sz in enumerate([17, 64, 128, 200]):
    S = set(random.sample(BOOLS, sz))
    a = {f: (1 if f in S else 0) for f in BOOLS}
    a[30163] = random.choice([0, 1]); a[11559] = random.choice([0, 1])
    configs.append(('rand|S|=%d' % sz, a))

for tag, extra in configs:
    seed = dict(C.BASE); seed.update(extra)
    v0 = E.forward(seed)
    ns = {'v': v0, '__builtins__': {}}
    cur = eval(H.acodes[TARGET], ns)
    offenders = []
    moved = 0
    for f in AFFINE:
        try:
            v1, _ = fast.apply_delta(v0, {f: v0[f] + 1})
        except Exception:
            continue
        ns1 = {'v': v1, '__builtins__': {}}
        try:
            new = eval(H.acodes[TARGET], ns1)
        except Exception:
            continue
        d = new - cur
        if d:
            moved += 1
            if d % P != 0:
                offenders.append((f, d))
    print("\n[%s] a20215 = %s... ; a20215 mod p = %s..." % (tag, str(cur)[:28], str(cur % P)[:28]), flush=True)
    print("    affine knobs that move a20215 at all: %d" % moved, flush=True)
    print("    of those, knobs whose step is NOT a multiple of p: %d" % len(offenders), flush=True)
    if offenders:
        print("    *** a20215 mod p IS MOVABLE -- p prime => 0 reachable. sec 3 REFUTED. ***", flush=True)
        for f, d in offenders[:6]:
            print("        x_%d step=%s (mod p = %s)" % (f, str(d)[:40], str(d % P)[:40]), flush=True)
    else:
        print("    every affine knob moves a20215 by an exact multiple of p", flush=True)
        print("    => a20215 mod p is invariant under ALL affine directions at this configuration", flush=True)
