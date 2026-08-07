"""Distribution of gap_p over the selector axis."""
import os, sys, json, glob
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REG = sys.argv[1] if len(sys.argv) > 1 else 'pinned'
pat = 'LATC_*.jsonl' if REG == 'consistent' else 'LAT_*.jsonl'
print('=== REGIME: %s   (%s) ===' % (REG, pat))
rows = []
for f in sorted(glob.glob(os.path.join(HERE, 'runs', pat))):
    for ln in open(f):
        try:
            rows.append(json.loads(ln))
        except Exception:
            pass
seen = {}
for r in rows:
    seen[r['tag']] = r
rows = list(seen.values())
ok = [r for r in rows if 'lat_gap_p' in r]
sk = [r for r in rows if 'lat_gap_p' not in r]

print('configurations measured: %d   (skipped/failed: %d)' % (len(ok), len(sk)))
for r in sk:
    print('   SKIP %-26s %s' % (r['tag'], r.get('note', '')))

print()
hdr = ('%-26s %-4s %-6s %-5s %-5s | %-4s %-4s | %-5s %-4s %-4s %-4s %-4s | %-4s %-4s | %s' %
       ('tag', 'live', 'score', '|R|', 'knobs', 'arQ', 'agp',
        'dim', 'lrQ', 'lrp', 'lgQ', 'lgp', 'rkq', 'gq', 'score<='))
print(hdr)
print('-' * len(hdr))
for r in sorted(ok, key=lambda r: (r['nlive'], r['R'])):
    print('%-26s %-4d %-6d %-5d %-5d | %-4d %-4d | %-5d %-4d %-4d %-4d %-4d | %-4s %-4s | %s' %
          (r['tag'], r['nlive'], r['score'], r['R'], r['knobs'],
           r['amb_rk_Q'], r['amb_gap_p'],
           r['lat_dim'], r['lat_rk_Q'], r['lat_rk_p'], r['lat_gap_Q'], r['lat_gap_p'],
           r.get('lat_rk_q'), r.get('lat_gap_q_ctl'), r.get('score_ub_p')))

print()
print('NOTE ON WHAT gap MEASURES.  [M|b] has exactly one column more than M, so')
print('    rk([M|b]) <= rk(M) + 1   and   gap in {0, 1}  BY CONSTRUCTION.')
print('gap = 1 is therefore the BOOLEAN "inconsistent"; it can never "widen", and a claim that')
print('it "stays exactly 1" is a claim of invariant INCONSISTENCY, not of an invariant magnitude.')
print('The quantitative p-statement is the DEFICIENCY rk_Q - rk_p, cross-checked against a')
print('control prime q of the same size: rk_Q - rk_q = 0 means the deficiency belongs to p.')
print()
print('=== DISTRIBUTION of lattice gap_p ===')
c = Counter(r['lat_gap_p'] for r in ok)
for k in sorted(c):
    print('  gap_p = %-3d : %3d configurations' % (k, c[k]))
print('=== DISTRIBUTION of lattice gap_Q ===')
c = Counter(r['lat_gap_Q'] for r in ok)
for k in sorted(c):
    print('  gap_Q = %-3d : %3d configurations' % (k, c[k]))
print('=== DISTRIBUTION of ambient gap_p (all knobs free) ===')
c = Counter(r['amb_gap_p'] for r in ok)
for k in sorted(c):
    print('  gap_p = %-3d : %3d configurations' % (k, c[k]))
print('=== mod-p rank DEFICIENCY rk_Q - rk_p on the lattice ===')
c = Counter(r['lat_deficiency'] for r in ok)
for k in sorted(c):
    print('  deficiency = %-3d : %3d configurations' % (k, c[k]))
print('=== CONTROL prime q (same size, unrelated): deficiency rk_Q - rk_q, and gap_q ===')
c = Counter(r.get('lat_deficiency_ctl') for r in ok)
for k in sorted(c, key=lambda x: (x is None, x)):
    print('  deficiency = %-4s : %3d configurations' % (k, c[k]))
c = Counter(r.get('lat_gap_q_ctl') for r in ok)
for k in sorted(c, key=lambda x: (x is None, x)):
    print('  gap_q      = %-4s : %3d configurations' % (k, c[k]))

print()
print('=== CROSS-TAB  gap_Q x gap_p  (the only p-SPECIFIC cell is gap_Q = 0, gap_p > 0) ===')
ct = Counter((r['lat_gap_Q'], r['lat_gap_p']) for r in ok)
for k in sorted(ct):
    mark = ''
    if k[0] == 0 and k[1] > 0:
        mark = '   <- genuine mod-p obstruction (solvable over Q, blocked mod p)'
    if k[0] == 0 and k[1] == 0:
        mark = '   <- TARGET: no obstruction at either level'
    if k[0] > 0:
        mark = '   <- already inconsistent over Q; gap_p carries no extra information'
    print('  gap_Q = %-3d gap_p = %-3d : %3d%s' % (k[0], k[1], ct[k], mark))
print('=== CROSS-TAB  gap_Q x gap_q(control prime) ===')
ct = Counter((r['lat_gap_Q'], r.get('lat_gap_q_ctl')) for r in ok)
for k in sorted(ct, key=lambda x: (x[0], -1 if x[1] is None else x[1])):
    print('  gap_Q = %-3d gap_q = %-4s : %3d' % (k[0], k[1], ct[k]))

z = [r for r in ok if r['lat_gap_p'] == 0]
print()
print('configurations with lattice gap_p = 0 (CONSISTENT mod p):', len(z))
for r in z:
    print('   %-26s score=%d |R|=%d dim=%d rkQ=%d rkp=%d gapQ=%d' %
          (r['tag'], r['score'], r['R'], r['lat_dim'], r['lat_rk_Q'], r['lat_rk_p'],
           r['lat_gap_Q']))
zq = [r for r in ok if r['lat_gap_Q'] == 0]
print('configurations with lattice gap_Q = 0 (solvable over Q):', len(zq),
      [r['tag'] for r in zq])
print()
print('score upper bound from the row-level p-obstruction (39033 - unzeroable_p):')
c = Counter(r.get('score_ub_p') for r in ok)
for k in sorted(c, key=lambda x: -(x or 0)):
    print('   <= %-6s : %3d configurations' % (k, c[k]))
print('degree overflow (a knob response of degree > %d anywhere): %d configurations' %
      (6, sum(1 for r in ok if r.get('degree_overflow'))))
