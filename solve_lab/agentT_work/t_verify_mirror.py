#!/usr/bin/env python3
"""AUDIT T37 -- verify the MIRROR by NUMBERS, against L's and F's published census.

Why this file exists, separately from t_rebuild*.sh: those scripts wrote `=== REBUILD DONE` at the
end and every stage was `python3 X.py | tail -3`.  A pipeline's exit status is its LAST command's,
so `tail` succeeded even when a stage crashed, `set -e` never fired, and the success marker could
be printed over a broken build.  (`set -o pipefail` is now set in all three -- but a status marker
must never be the evidence.)  This script is the evidence: every quantity is asserted against a
closed form and the process EXITS NONZERO on any mismatch.  It prints no success word it has not
earned.

  python3 t_verify_mirror.py            structure only, ~10 s
  python3 t_verify_mirror.py --deep     + reload close_T8.json through L's engine, ~90 s
"""
import os, sys, pickle, collections

T = os.path.dirname(os.path.abspath(__file__))
MF = os.path.join(T, 'mirror', 'F')
ML = os.path.join(T, 'mirror', 'L')
fails = []

CHECKS = []

def chk(name, got, want):
    CHECKS.append(name)
    ok = (got == want)
    print('   %-46s %-30s %s' % (name, str(got)[:30], 'OK' if ok else 'MISMATCH, want %s' % (want,)))
    if not ok:
        fails.append((name, got, want))

print('== F chain ==')
d = pickle.load(open(os.path.join(MF, 'circ4.pkl'), 'rb'))
chk('circ4: atoms', len(d['atoms']), 39033)
chk('circ4: equations', len(d['eqrows']), 39033)
s = pickle.load(open(os.path.join(MF, 'sched.pkl'), 'rb'))
chk('sched: scheduled defs', len(s['order']), 30001)

print('== L chain ==')
ors = pickle.load(open(os.path.join(ML, 'ors.pkl'), 'rb'))
chk('global: OR nodes', len(ors), 383)
ot = pickle.load(open(os.path.join(ML, 'ortree2.pkl'), 'rb'))
# `tree` holds internal nodes AND leaves; the published census is "OR nodes 254, leaves 256",
# so the closed form is the DECOMPOSITION, not the total.  (My first version of this file
# asserted 254 against len(tree) and flagged a mismatch that was my own mis-keying -- rule 7:
# separate "this number is wrong" from "this result is wrong".)
chk('ortree2: internal OR nodes', sum(1 for v in ot['tree'].values() if v), 254)
chk('ortree2: leaves', sum(1 for v in ot['tree'].values() if not v), 256)
chk('ortree2: tree total = 254+256', len(ot['tree']), 510)
chk('ortree2: one selector per OR node', len(ot['selmap']), 254)
H = pickle.load(open(os.path.join(ML, 'handles.pkl'), 'rb'))
# handles.py PRINTS 3,707 -- that is the probe count (every free var whose atom-deltas are all
# divisible by p).  The stored `handle` map is L's 3,681 census.  Both are real and they
# reconcile as 3,707 - 33 + 7 = 3,681 (RESUME_T sixth pass).  Assert both, and the partition.
chk('handles: handle free vars (L census)', len(H['handle']), 3681)
chk('handles: handle + value == free vars', len(set(H['handle']) | set(H['value'])), 8747)
M = pickle.load(open(os.path.join(ML, 'full_model.pkl'), 'rb'))
chk('buildall: nodes', len(M['NODE']), 383)
chk('buildall: live leaves', len(M['live']), 256)
chk('buildall: dead leaves', len(M['dead']), 128)
chk('buildall: pins', len(M['PIN']), 256)
C = pickle.load(open(os.path.join(ML, 'calib2.pkl'), 'rb'))
oh = collections.Counter(C['ORIENT'].values())
chk('calib2: orient hist', sorted(oh.items(), key=lambda kv: str(kv[0])),
    [(0, 67), (1, 188), ('DEAD', 128)])
SL = pickle.load(open(os.path.join(ML, 'slopes.pkl'), 'rb'))
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
chk('slopes: handle atoms', len(SL), 3681)
chk('slopes: slopes not divisible by p', sum(1 for v in SL.values() if v and v % p), 0)
hist = collections.Counter((abs(v)//p if v else 0) for v in SL.values())
chk('slopes: c == 1', hist[1], 2747)
chk('slopes: zero slope', hist[0], 7)
chk('slopes: c > 1  (THE 927)', sum(n for c, n in hist.items() if c > 1), 927)

if '--deep' in sys.argv:
    print('== end-to-end: reload a KNOWN state through L\'s engine ==')
    sys.path.insert(0, T)
    import t_close2w as CC
    import json
    vv = [0]*CC.NV
    for k, val in json.load(open(os.path.join(T, 'close_T8.json'))).items():
        vv[int(k[2:])] = int(val)
    CC.relift(vv)
    r = CC.E.run(vv)
    chk('close_T8.json nonzero atoms', sum(1 for x in r if x), 3)
    chk('L engine residual atoms', len(CC.E.res), 9032)

print()
if fails:
    print('MIRROR IS NOT FAITHFUL -- %d mismatch(es): %s' % (len(fails), [f[0] for f in fails]))
    sys.exit(1)
print('mirror verified against %d closed-form quantities, 0 mismatches' %
      (len(CHECKS)))
