"""Full-coverage sweep over detach sets, using the measured collapse of the placement space.

Measured: across placements the region R takes only two values, and the zero-collateral knob
list and the knob-response matrix M depend only on R -- only the row-target vector b varies,
and b takes very few distinct values.  So price a placement by (R, b) and cache; the expensive
part (probing every knob to build M) then runs once per NEW signature instead of once per set.

A random audit re-does the full build on a fraction of cache HITS and checks that the knob list
and M really do match the cached ones.  The audit rate and result are reported, because the
speed-up is only sound if that holds.

usage: python3 sweepF.py <out.jsonl> <k> [shard nshard] [audit_rate]
"""
import sys, json, time, os, itertools, hashlib, random
import optN
from optN import make, build, POOL, inner, atom_eqs
import ev
import zsolve

out = sys.argv[1]
k = int(sys.argv[2])
shard = int(sys.argv[3]) if len(sys.argv) > 3 else 0
nsh = int(sys.argv[4]) if len(sys.argv) > 4 else 1
audit_rate = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
random.seed(1234 + shard)

cache = {}          # (R, bhash) -> record
mhash_by_R = {}     # R -> hash of (knobs, M)
audits = [0, 0]     # done, mismatched


def region_of(st):
    NZ = set(st.nz())
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    return tuple(sorted(R))


def price_fast(D):
    st = make(list(D))
    R = region_of(st)
    b = tuple(inner(st, e) for e in R)
    bh = hashlib.sha256(repr(b).encode()).hexdigest()[:16]
    key = (R, bh)
    hit = key in cache
    do_audit = hit and random.random() < audit_rate
    if hit and not do_audit:
        return dict(cache[key]), 'hit'
    d = build(st)
    assert tuple(d['R']) == R and tuple(d['b']) == b
    mh = hashlib.sha256(repr((d['knobs'], d['M'])).encode()).hexdigest()[:16]
    if R in mhash_by_R and mhash_by_R[R] != mh:
        print('!!! M DIFFERS for the same region at D=%s' % (list(D),), flush=True)
    mhash_by_R.setdefault(R, mh)
    outside = d['outside']
    opt, rows, exh, tests = zsolve.max_zero_rows(d['M'], d['b'], d['n'], len(R))
    fail = len(R) - opt + outside
    rec = dict(R=len(R), knobs=d['n'], opt=opt, outside=outside, failing=fail,
               score=39033 - fail, exh=str(exh), mh=mh, bh=bh)
    if do_audit:
        audits[0] += 1
        old = cache[key]
        same = (old['opt'] == rec['opt'] and old['failing'] == rec['failing']
                and old['mh'] == rec['mh'])
        if not same:
            audits[1] += 1
            print('!!! AUDIT MISMATCH at D=%s cached=%s fresh=%s' % (list(D), old, rec), flush=True)
        return dict(rec), 'audit'
    cache[key] = rec
    return dict(rec), 'miss'


sets = [c for i, c in enumerate(itertools.combinations(POOL, k)) if i % nsh == shard]
done = set()
if os.path.exists(out):
    for ln in open(out):
        try:
            done.add(tuple(json.loads(ln)['D']))
        except Exception:
            pass
todo = [s for s in sets if s not in done]
print('k=%d sets=%d done=%d todo=%d audit_rate=%.3f' % (k, len(sets), len(done), len(todo),
                                                        audit_rate), flush=True)
f = open(out, 'a')
t0 = time.time()
best = 0
nmiss = 0
for i, D in enumerate(todo):
    rec, how = price_fast(D)
    if how == 'miss':
        nmiss += 1
    rec['D'] = list(D)
    rec['how'] = how
    f.write(json.dumps(rec) + '\n')
    if rec['score'] > best:
        best = rec['score']
    if rec['score'] > 39026:
        print('*** BEATS 39026: %s score=%d ***' % (list(D), rec['score']), flush=True)
        f.flush()
    if (i + 1) % 2000 == 0:
        f.flush()
        el = time.time() - t0
        print('%d/%d  %.4fs/set  eta %.1f min  signatures=%d audits=%d/%d bad  best=%d'
              % (i + 1, len(todo), el / (i + 1), (len(todo) - i - 1) * el / (i + 1) / 60,
                 len(cache), audits[0], audits[1], best), flush=True)
f.flush()
f.close()
print('DONE k=%d shard=%d: %d sets, %d distinct (R,b) signatures, %d audits with %d mismatches, '
      'best=%d, %.1f min'
      % (k, shard, len(todo), len(cache), audits[0], audits[1], best, (time.time() - t0) / 60),
      flush=True)
