"""AMBIENT p-deficiency of the region response, across the structural selector axis.

The lattice measurement (pselsum.py) answers "is the region BLOCKED mod p".  This answers the
prior question the task actually names — *is the region response rank-deficient mod p* — over the
whole knob space, where it is cheap enough to run on every configuration.

Two quantities, and the control prime is what makes them mean anything:
    deficiency_p   = rk_Q(M) - rk_p(M)
    deficiency_q   = rk_Q(M) - rk_q(M)   with q a prime of the SAME SIZE, unrelated to the frame
A generic prime loses no rank, so deficiency_q = 0 is the null result; deficiency_p > 0 with
deficiency_q = 0 is a fact about p.

    pq_knobs_p / knobs_live = the fraction of knobs that actually move the region whose ENTIRE
    region column is 0 mod p — knobs that move the region only in multiples of p.  That is the
    MECHANISM behind the deficiency, not a restatement of it.
"""
import os, sys, json, glob
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
rows = {}
for f in sorted(glob.glob(os.path.join(HERE, 'runs', sys.argv[1] if len(sys.argv) > 1
                                       else 'pselrank_amb.jsonl'))):
    for ln in open(f):
        try:
            r = json.loads(ln)
            rows[r['tag']] = r
        except Exception:
            pass
ok = [r for r in rows.values() if 'amb_rk_q' in r]
sk = [r for r in rows.values() if 'amb_rk_q' not in r]
print('configurations: %d measured, %d skipped' % (len(ok), len(sk)))
for r in sk:
    print('   SKIP %-26s %s' % (r['tag'], r.get('note', '')))
print()
h = ('%-26s %-4s %-6s %-5s %-5s | %-5s %-5s %-5s | %-6s %-6s | %-5s %-5s %s' %
     ('tag', 'live', 'score', '|R|', 'knobs', 'rk_Q', 'rk_p', 'rk_q',
      'def_p', 'def_q', 'pq_p', 'pq_q', 'live_cols'))
print(h)
print('-' * len(h))
for r in sorted(ok, key=lambda r: (r['nlive'], r['R'])):
    print('%-26s %-4d %-6d %-5d %-5d | %-5d %-5d %-5d | %-6d %-6d | %-5d %-5d %d' %
          (r['tag'], r['nlive'], r['score'], r['R'], r['knobs'],
           r['amb_rk_Q'], r['amb_rk_p'], r['amb_rk_q'],
           r['amb_deficiency'], r['amb_deficiency_ctl'],
           r['pq_knobs_p'], r['pq_knobs_ctl'], r['knobs_live']))
print()
print('=== DISTRIBUTION: deficiency mod p  (rk_Q - rk_p) ===')
c = Counter(r['amb_deficiency'] for r in ok)
print('   zero deficiency (NOT rank-deficient mod p): %d of %d' % (c.get(0, 0), len(ok)))
print('   min %d   max %d   median %d' %
      (min(c), max(c), sorted(r['amb_deficiency'] for r in ok)[len(ok) // 2]))
print('=== DISTRIBUTION: deficiency mod the CONTROL prime q (rk_Q - rk_q) ===')
c = Counter(r['amb_deficiency_ctl'] for r in ok)
for k in sorted(c):
    print('   deficiency_q = %-4d : %3d configurations' % (k, c[k]))
print('=== p-QUANTISED KNOBS: whole region column = 0 mod p, of the knobs that move R ===')
fr = [(r['pq_knobs_p'], r['knobs_live'], r['tag']) for r in ok]
tot_p = sum(a for a, b, t in fr)
tot_l = sum(b for a, b, t in fr)
print('   mod p : %d of %d live knob columns (%.1f%%)' % (tot_p, tot_l, 100.0 * tot_p / tot_l))
print('   mod q : %d of %d live knob columns (control)' %
      (sum(r['pq_knobs_ctl'] for r in ok), tot_l))
print('   per-configuration fraction: min %.2f  max %.2f' %
      (min(a / b for a, b, t in fr if b), max(a / b for a, b, t in fr if b)))
print('=== ambient gap (is b in the column span?) ===')
print('   gap_Q distribution:', dict(Counter(r['amb_gap_Q'] for r in ok)))
print('   gap_p distribution:', dict(Counter(r['amb_gap_p'] for r in ok)))
print('   gap_q distribution:', dict(Counter(r['amb_gap_q'] for r in ok)))
