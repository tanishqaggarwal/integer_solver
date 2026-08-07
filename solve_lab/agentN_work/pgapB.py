"""Part B of the rank-gap experiment, bounded: verified assignments as base configurations.

Part A (the 16 detach states) is in pgap.py and is the definitive statement.  Here every verified
assignment on disk is loaded as a base state and priced the same way, with an explicit size guard so
one huge region cannot stall the sweep — anything skipped is reported as skipped, not omitted.
"""
import os, sys, json, glob, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import pgap
from pgap import price, show, HDR, fr
from frameB import State
from optN import atom_eqs

MAXR, MAXK = 40, 220
files = sorted(glob.glob(os.path.join(HERE, '*.json')) +
               glob.glob(os.path.join(HERE, '..', 'best', '*.json')))
print(HDR, flush=True)
out = []
for f in files:
    try:
        W = json.load(open(f))
    except Exception:
        continue
    if not isinstance(W, dict) or len(W) < 1000:
        continue
    tag = os.path.basename(f)[:32]
    try:
        v = [0] * 38748
        for kk, val in W.items():
            v[int(kk[2:]) if str(kk).startswith('x_') else int(kk)] = int(val)
        fv = {u: v[u] for u in fr.free if v[u] != 0}
        st = State(fr, fv)
    except Exception as e:
        print('%-34s SKIP (load: %s)' % (tag, str(e)[:40]), flush=True); continue
    NZ = set(st.nz())
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    if len(R) > MAXR:
        print('%-34s SKIP (region |R|=%d > %d, frame score %d)' % (tag, len(R), MAXR, st.score()),
              flush=True)
        out.append(dict(tag=tag, skipped='|R|=%d' % len(R), frame_score=st.score()))
        continue
    t0 = time.time()
    try:
        r = price(st, tag)
    except Exception as e:
        print('%-34s SKIP (price: %s)' % (tag, str(e)[:40]), flush=True); continue
    r['frame_score'] = st.score(); r['secs'] = round(time.time() - t0, 1)
    show(r); out.append(r)

pr = [r for r in out if 'gap_p' in r]
print('\npriced %d, skipped %d' % (len(pr), len(out) - len(pr)), flush=True)
print('distinct mod-p gaps : %s' % sorted(set(r['gap_p'] for r in pr)), flush=True)
print('distinct over-Q gaps: %s' % sorted(set(r['gap_Q'] for r in pr)), flush=True)
if pr:
    b = max(pr, key=lambda r: r['score'])
    print('best score priced   : %d (%s)' % (b['score'], b['tag']), flush=True)
json.dump(out, open(os.path.join(HERE, 'runs', 'pgapB.json'), 'w'), indent=1)
