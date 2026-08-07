#!/usr/bin/env python3
"""agent AE -- assemble the family results table from the run artefacts.

Every row's evidence is the engine's own DONE line: a jump count, a DP count, and the
ratio of the DP count to the closed form jumps/2^dpbits.  A row without a DONE line is
reported as NOT RUN, never as a miss.
"""
import json, os, math, re, glob

N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
LGN = 256.0

def conf(jumps, R):
    """P(a hit inside the interval would have been found), exponential model, mean 2*sqrt(L)"""
    return 1.0 - math.exp(-jumps / (2.0 * 2 ** (R / 2.0)))

def parse_done(path):
    if not os.path.exists(path): return None
    txt = open(path).read()
    m = None
    for line in txt.splitlines():
        if line.startswith('DONE'): m = line
    if m is None:
        st = [l for l in txt.splitlines() if l.startswith('STATUS')]
        return dict(partial=True, line=st[-1] if st else None,
                    jumps=int(re.search(r'jumps=(\d+)', st[-1]).group(1)) if st else 0)
    d = dict(partial=False, line=m)
    for kv in m.split()[1:]:
        k, v = kv.split('=')
        try: d[k] = float(v) if ('.' in v) else int(v)
        except ValueError: d[k] = v
    return d

rows = []

# --- headline
for tag, R, f in (('R=58 magnitude k0 < 2^58', 58, 'head58.err'),
                  ('R=60 magnitude (abandoned)', 60, 'head60_PARTIAL_ABANDONED.txt'),
                  ('R=64 magnitude (abandoned)', 64, 'head64_PARTIAL_ABANDONED.txt')):
    d = parse_done(f)
    if not d: continue
    j = d.get('jumps', 0)
    rows.append(dict(family=tag, R=R, n=1, jumps=j, size=2.0 ** R,
                     dpcf=(d['dps'] / (d['jumps'] / 2.0 ** 14)) if (not d['partial'] and 'dps' in d) else None,
                     complete=not d['partial'], conf=conf(j, R), hit=False))

# --- tiers
for tier, fn in (('const', 'res_const.json'), ('orbit', 'res_orbit.json'), ('window', 'res_window.json')):
    if not os.path.exists(fn): continue
    rs = json.load(open(fn))
    done = [r for r in rs if r.get('jumps') and r['rc'] == 0]
    R = rs[0]['R']
    tot = sum(r['jumps'] for r in done)
    hits = [r for r in rs if r.get('hit')]
    cf = [r['dp_over_closedform'] for r in done if r.get('dp_over_closedform')]
    rows.append(dict(family='%s tier (%d families, R=%d)' % (tier, len(rs), R), R=R, n=len(done),
                     jumps=tot, size=len(done) * 2.0 ** R,
                     dpcf=(min(cf), max(cf)) if cf else None,
                     complete=(len(done) == len(rs)),
                     conf=conf(min(r['jumps'] for r in done), R) if done else 0.0,
                     hit=bool(hits)))

# --- quotient
for f in sorted(glob.glob('res_quot_*.json')):
    q = json.load(open(f))
    m = int(f.split('_')[-1].split('.')[0])
    sz = 6 * (6 / math.pi ** 2) * (2.0 ** m) ** 2
    rows.append(dict(family='quotient a*b^-1, m=%d (EXHAUSTIVE)' % m, R=None, n=1,
                     jumps=2 * 2 ** m, size=sz, dpcf='n/a', complete=(q['rc'] == 0),
                     conf=1.0, hit=bool(q['verified'])))

print('| family | keys covered | prior P(k0 in F) | ops | exclusion confidence | outcome |')
print('|---|---|---|---|---|---|')
tot = 0.0
for r in rows:
    if 'abandoned' in r['family']:
        st = 'NOT RUN (partial, %.0f%% conf) -- reported as not-run' % (100 * r['conf'])
    elif not r['complete']:
        st = 'INCOMPLETE'
    elif r['hit']:
        st = '*** HIT ***'
    else:
        st = 'exhausted-no-hit' if 'EXHAUSTIVE' in r['family'] else 'miss'
        tot += r['size']
    print('| %s | 2^%.2f | 2^%.1f | 2^%.1f | %s | %s |' % (
        r['family'], math.log2(r['size']), math.log2(r['size']) - LGN,
        math.log2(max(r['jumps'], 1)),
        ('%.1f%%' % (100 * r['conf'])) if r['conf'] < 1 else '100% (exhaustive)', st))
print()
print('TOTAL keys excluded by rows that completed: 2^%.2f  (prior mass 2^%.1f)'
      % (math.log2(tot) if tot else 0, (math.log2(tot) - LGN) if tot else 0))
print('For comparison: agent X weight<=9 covered 2^53.38; agent Y complement 2^53.38.')
