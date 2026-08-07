#!/usr/bin/env python3
"""Per-offset results table: exact candidate counts, exhaustion status, hits."""
import json, os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'aa_offsets.json')))
man = {x['tag']: x for x in D['manifest']}
pre = sys.argv[1] if len(sys.argv) > 1 else 'd'
EXP = {1: 512, 2: 130560, 3: 22108160, 4: 2796682240}
rows = []
for tag, x in man.items():
    f = os.path.join(HERE, 'runs', 'r_%s_%s.txt' % (pre, tag))
    done, n, z, t, hits = {}, 0, 0, 0.0, 0
    if os.path.exists(f):
        for ln in open(f):
            p = ln.split()
            if p and p[0] == 'HIT': hits += 1
            elif p and p[0] == 'DONE':
                kv = dict(q.split('=', 1) for q in p if '=' in q)
                sz = int(kv['sz']); done[sz] = int(kv['n']); z += int(kv['zero'])
                t += float(p[-1].rstrip('s'))
    rows.append((x['tier'], -(x['reach'] or 0), tag, x['reach'], done, z, t, hits))
rows.sort()
print('%-11s %-4s %-5s %-30s %-6s %-6s %-8s %s' %
      ('tag', 'tier', 'reach', 'scan sizes complete (exact n)', 'm<=', 'zero', 'sec', 'HITS'))
tot_c = tot_t = tot_h = 0
for tier, _, tag, reach, done, z, t, hits in rows:
    ok = [b for b in sorted(done) if done[b] == EXP[b]]
    bad = [b for b in sorted(done) if done[b] != EXP[b]]
    mmax = 4 + max(ok) if ok and set(ok) >= set(range(1, max(ok) + 1)) else 0
    tot_c += sum(done.values()); tot_t += t; tot_h += hits
    print('%-11s %-4d %-5s %-30s %-6s %-6d %-8.1f %d%s' %
          (tag, tier, reach, ','.join('b%d:%d' % (b, done[b]) for b in sorted(done)) or '-',
           mmax or '-', z, t, hits, '  MISMATCH %s' % bad if bad else ''))
print('\ntotals: %d scan candidates, %.1f s of engine time, %d hits, %d offsets'
      % (tot_c, tot_t, tot_h, len(rows)))
print('expected false positives over this sweep: %.3f'
      % (tot_c * 1409460736 / 2.0**64))
