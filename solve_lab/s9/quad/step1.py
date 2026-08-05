"""Step 1: catalogue the residual on each quadrant branch."""
import sys, time, collections
sys.path.insert(0, 'quad')
from common import *

t0 = time.time()
v0 = H.load_assignment(BEST)
CODES, _ = H.load_equations()

print('=== BASELINE (best partial) ===')
print('x_2081 =', v0[2081], ' x_24601 =', v0[24601], ' x_15298 =', v0[15298],
      ' x_7715 =', v0[7715], ' x_34554 =', v0[34554], ' x_4287 =', v0[4287])
nz0 = nz_checks(v0)
print('nonzero residual checks:', len(nz0), nz0)
print('nonzero raw atoms:', [a for a in range(len(polys)) if evalpoly(polys[a], v0) != 0])
f0 = H.evaluate(CODES, v0)
print(f'FULL harness: {len(CODES)-len(f0)}/{len(CODES)}  ({len(f0)} failing) {f0}')
# validate the atom2eq shortcut
live = [a for a in range(len(polys)) if evalpoly(polys[a], v0) != 0]
cand = eqs_of(live)
ff = H.evaluate(CODES, v0, cand)
print(f'SHORTCUT: {len(ff)} failing, matches full = {sorted(ff)==sorted(f0)}')
print(f'({time.time()-t0:.0f}s)')


def classify(a, v):
    s = src[a]
    r = evalpoly(resid_poly[a], v)
    vs = sorted(set(u for m in resid_poly[a] for u in m))
    tag = 'SQUARE-ROOT' if a in roots else 'plain'
    return (a, tag, r, len(vs), s)


def report(label, seeds):
    print('\n' + '=' * 78)
    print('BRANCH', label, 'seeds', seeds)
    v = list(v0)
    ch, steps = ripple(v, seeds)
    print(f'  ripple changed {len(ch)} vars in {steps} steps')
    print(f'  x_15298={v[15298]} x_7715={v[7715]} x_34554={v[34554]} x_19247={v[19247]}')
    ng = nz_gates(v)
    print(f'  nonzero GATE atoms (should be 0): {len(ng)} {ng[:10]}')
    nz = nz_checks(v)
    print(f'  nonzero residual checks: {len(nz)}')
    for a in nz:
        aid, tag, r, nv, s = classify(a, v)
        print(f'   atom {aid:6d} [{tag:11s}] |resid|={r.bit_length():4d}b  rmodp={"0" if r % P == 0 else "nz"}  nvars={nv}')
        print(f'        src: {s[:200]}')
    liv = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    cand = eqs_of(liv)
    ff = H.evaluate(CODES, v, cand)
    print(f'  live raw atoms {len(liv)} -> {len(cand)} candidate eqs -> FAILING {len(ff)}')
    return v, nz, ff


for lbl, sd in [('x_2081=0', {2081: 0}), ('x_24601=0', {24601: 0}), ('both=0', {2081: 0, 24601: 0})]:
    report(lbl, sd)
print(f'\ntotal {time.time()-t0:.0f}s')
