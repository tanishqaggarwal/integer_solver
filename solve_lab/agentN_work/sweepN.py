"""Resumable corrected-stage-B sweep over detach sets.

usage: python3 sweepN.py <out.jsonl> <k> <shard> <nshard>     # all k-subsets of POOL
       python3 sweepN.py <out.jsonl> list <file.json>          # explicit list of sets
Writes one JSON record per line; re-running skips sets already present.
"""
import sys, json, time, os, itertools
import optN
from optN import POOL, price_D

out = sys.argv[1]
mode = sys.argv[2]

if mode == 'list':
    sets = [tuple(sorted(s)) for s in json.load(open(sys.argv[3]))]
else:
    k = int(mode)
    shard = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    nsh = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    sets = [c for i, c in enumerate(itertools.combinations(POOL, k)) if i % nsh == shard]

done = set()
if os.path.exists(out):
    for ln in open(out):
        try:
            done.add(tuple(json.loads(ln)['D']))
        except Exception:
            pass
todo = [s for s in sets if s not in done]
print('sets=%d done=%d todo=%d' % (len(sets), len(done), len(todo)), flush=True)

f = open(out, 'a')
t0 = time.time()
best = 0
for i, D in enumerate(todo):
    try:
        r = price_D(list(D), want_lin=False)
    except Exception as e:
        f.write(json.dumps({'D': list(D), 'err': str(e)[:200]}) + '\n')
        f.flush()
        continue
    rec = {'D': list(D), 'R': r['R'], 'S': r['S'], 'knobs': r['knobs'], 'rank': r['rank'],
           'z0': r['z0'], 'opt': r['opt'], 'outside': r['outside'], 'failing': r['failing'],
           'score': r['score'], 'exh': r['exhaustive'], 'lin': r['lin']}
    if r['score'] > 39026:
        rec['rows'] = r['rows']
        rec['knoblist'] = r['knoblist']
        print('*** BEATS 39026: %s score=%d ***' % (list(D), r['score']), flush=True)
    f.write(json.dumps(rec) + '\n')
    if r['score'] > best:
        best = r['score']
    if (i + 1) % 200 == 0:
        f.flush()
        el = time.time() - t0
        print('%d/%d  %.3fs/set  eta %.1f min  best=%d' %
              (i + 1, len(todo), el / (i + 1), (len(todo) - i - 1) * el / (i + 1) / 60, best),
              flush=True)
f.flush()
f.close()
print('DONE %d sets in %.1f min, best=%d' % (len(todo), (time.time() - t0) / 60, best), flush=True)
