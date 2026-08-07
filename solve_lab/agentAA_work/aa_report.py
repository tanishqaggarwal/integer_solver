#!/usr/bin/env python3
"""FINAL REPORT -- strictly count-based.

Per the coordinator's rule (from agent Y's `yorbit.status` incident): a status marker is a
claim, a count checked against a closed-form expectation is evidence.  This script reads ONLY
the engine's own `DONE ... n=<count>` lines and checks every count against C(256,b)*2^b.  The
`OFFSET_DONE` / `SHARDED_DONE` / `SHARD<s>` markers written by the shell drivers are IGNORED
here -- they were echoed without testing an exit code and are therefore worthless as evidence.
Anything without counts is printed as NEVER RUN, not omitted.
"""
import json, os, sys
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'aa_offsets.json')))
EXP = {b: comb(256, b) * 2**b for b in (1, 2, 3, 4)}

def counts(path):
    """-> {sz: [n, n, ...]} one entry per DONE line actually emitted by the engine"""
    d = {}
    if os.path.exists(path):
        for ln in open(path):
            p = ln.split()
            if p and p[0] == 'DONE':
                kv = dict(q.split('=', 1) for q in p if '=' in q)
                d.setdefault(int(kv['sz']), []).append(int(kv['n']))
    return d

def zeros(path):
    z = 0
    if os.path.exists(path):
        for ln in open(path):
            p = ln.split()
            if p and p[0] == 'DONE':
                kv = dict(q.split('=', 1) for q in p if '=' in q)
                z += int(kv['zero'])
    return z

def hits(path):
    return sum(1 for ln in open(path) if ln.startswith('HIT')) if os.path.exists(path) else 0

rows = []
for x in sorted(D['manifest'], key=lambda z: (z['tier'], -(z['reach'] or 0))):
    t = x['tag']
    f7 = os.path.join(HERE, 'runs', 'r_d_%s.txt' % t)          # full-table, m<=7
    fs = os.path.join(HERE, 'runs', 'rs_%s.txt' % t)           # 8 shard passes, m<=7
    f6 = os.path.join(HERE, 'runs6', 'r6_d_%s.txt' % t)        # a<=3 table, m<=6
    c7, cs, c6 = counts(f7), counts(fs), counts(f6)
    # m<=7 by the monolithic table: exactly one exhaustive DONE per b=1,2,3
    ok7 = all(c7.get(b, []).count(EXP[b]) >= 1 for b in (1, 2, 3))
    # m<=7 by shard passes: the 8 passes partition the table (shard = key>>61), so a size is
    # exhausted only when ALL EIGHT passes emitted the full candidate count for it
    oks = all(cs.get(b, []).count(EXP[b]) == 8 for b in (1, 2, 3))
    ok6 = all(c6.get(b, []).count(EXP[b]) >= 1 for b in (1, 2, 3))
    m = 7 if (ok7 or oks) else (6 if ok6 else 0)
    ev = ('full-table' if ok7 else ('8/8 shard passes' if oks else
          ('%d/8 shard passes' % min(cs.get(b, []).count(EXP[b]) for b in (1, 2, 3))
           if cs else ('a<=3 table' if ok6 else 'NEVER RUN'))))
    if m == 6 and not ok7 and not oks: ev = 'a<=3 table'
    h = hits(f7) + hits(fs) + hits(f6)
    z = zeros(f7) + zeros(fs) + zeros(f6)
    n = sum(sum(v) for v in list(c7.values()) + list(cs.values()) + list(c6.values()))
    rows.append((t, x['tier'], x['reach'], m, ev, n, z, h))

print('%-11s %-5s %-6s %-5s %-20s %-14s %-5s %s' %
      ('tag', 'tier', 'reach', 'm<=', 'evidence', 'candidates', 'zero', 'HITS'))
for r in rows:
    print('%-11s %-5d %-6s %-5s %-20s %-14d %-5d %d' %
          (r[0], r[1], r[2], (r[3] or 'NONE'), r[4], r[5], r[6], r[7]))

n7 = sum(1 for r in rows if r[3] >= 7); n6 = sum(1 for r in rows if r[3] >= 6)
nn = sum(1 for r in rows if r[3] == 0)
TOT = sum(r[5] for r in rows)
print('\nEXHAUSTED at m<=7 : %d of %d offsets' % (n7, len(rows)))
print('EXHAUSTED at m<=6 : %d of %d offsets  (all of the above included)' % (n6, len(rows)))
print('NEVER RUN         : %d' % nn)
print('total scan candidates (counted, not claimed): %d' % TOT)
print('degenerate dx=0 events: %d' % sum(r[6] for r in rows))
print('HITS: %d      expected false positives: %.3f' % (sum(r[7] for r in rows),
      TOT * 1409460736 / 2.0**64))
