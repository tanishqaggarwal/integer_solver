"""S11 step 93: boolean-flip census FROM THE ADVICE-SOLVED STATE.

§141 leaves the branch as the only thing in the instance that is not pinned to its own
constants.  The lab has censused boolean flips before and priced them at >= 8 against a
gain of 0 -- but that was measured at states where the advice congruences were still
unsolved and the residual was a different object.  Redo it from the 39,013 attractor,
where every advice value sits at its pin and the residual is exactly A != 0, B != 0.

For every free input that is 0/1 at the base state, flip it, forward-evaluate, and
record the score, the surviving nonzero checks, and what happens to the selector
x15298 and to A and B.  Chunked and resumable.

Usage: boolcensus.py START END [state.json]
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from chunk import sweep, load
P = ad.P
src = sys.argv[3] if len(sys.argv) > 3 else 'PIN_39013.json'
tag = 'bool_' + os.path.basename(src).replace('.json', '')
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
FREE = [t for t in range(L.NVARS) if t not in L.definer]
CAND = [t for t in FREE if base[t] in (0, 1)]
print('%s: score %d; %d boolean-valued free inputs' % (src, BASE, len(CAND)),
      flush=True)
print('base: x15298 = %d, A = %d, B = %d'
      % (base[15298], base[35389] % P, base[6671] % P), flush=True)


def evaluate(u):
    out = {'u': u, 'score': -1}
    for val in (0, 1):
        if val == base[u]:
            continue
        v = list(base)
        v[u] = val
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(av))
        nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
        if s > out['score']:
            out.update(score=s, val=val, nz=len(nz), sel=v[15298],
                       A=int(v[35389] % P == 0), B=int(v[6671] % P == 0))
        if s > BASE:
            T.save(v, os.path.join(HERE, 'BC_%d_x%d.json' % (s, u)))
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(CAND)
sweep(tag, CAND, evaluate, start, min(end, len(CAND)), keyfn=str, budget=540)
rs = load(tag)
if rs:
    rs.sort(key=lambda r: -r['score'])
    print('\ntop flips from %d:' % BASE)
    for r in rs[:15]:
        print('   x%-6d -> %-3s score %-6d checks %-4s selector %-3s A=0:%s B=0:%s'
              % (r['u'], r.get('val'), r['score'], r.get('nz'), r.get('sel'),
                 r.get('A'), r.get('B')))
    sel0 = [r for r in rs if r.get('sel') == 0]
    print('\nflips that drive the selector x15298 to 0: %d' % len(sel0))
    for r in sel0[:10]:
        print('   x%-6d score %d' % (r['u'], r['score']))
