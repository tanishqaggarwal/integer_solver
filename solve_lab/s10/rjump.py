"""S11 step 75: RESIDUE JUMPS -- a move class none of the vetoes ever priced.

Two facts settled this session change what a "move" is:

  * fidelity.py: no gate ever breaks under a large move (every pivot is +-1), so the
    47% misprediction rate is pure SECOND-ORDER content.  The map free-inputs ->
    checks is an honest polynomial, and every linear veto in this lab priced only
    its tangent space.
  * gfix.py: a failing check c that depends on a free input u can be zeroed mod p
    EXACTLY by setting u <- u + delta with delta = -c/(dc/du) mod p, and the handle
    then absorbs the quotient over Z.  That is not a step in the tangent space; it
    is a jump of a full residue class, and its true effect is only visible by
    evaluating.

So: enumerate every (failing check, free input in its support) pair, take the exact
residue jump, forward-evaluate, and score.  Nothing is predicted -- everything is
measured.  Chunked and resumable.

Usage: rjump.py START END [state.json]
"""
import os, sys, time, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
import suppfree
from chunk import sweep, load
P = ad.P
src = sys.argv[3] if len(sys.argv) > 3 else 'mod9118_0.json'
tag = 'rjump_' + os.path.basename(src).replace('.json', '')

base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
BAD = [a for a in range(L.NA) if a not in L.atom_out and bav[a]]
FREE = set(t for t in range(L.NVARS) if t not in L.definer)
print(f'{src}: score {BASE}; failing checks {BAD}', flush=True)

idx, freelist, vs = suppfree.build(bm)
CAND = []
for c in BAD:
    m = suppfree.atom_supp(c, bm, vs)
    sup = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]
    col = {}
    for u in sup:
        col[u] = jac_column(u, base, bm, [c]).get(c, 0)
    for u in sup:
        if col[u]:
            CAND.append((c, u))
print(f'{len(CAND)} (check, free input) residue jumps available', flush=True)


def evaluate(spec):
    c, u = spec
    g = jac_column(u, base, bm, [c]).get(c, 0)
    if not g:
        return {'score': -1, 'c': c, 'u': u}
    delta = (-bav[c] % P) * pow(g, -1, P) % P
    out = {'c': c, 'u': u, 'score': -1}
    for d in (delta, delta - P):
        v = list(base)
        v[u] = v[u] + d
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(av))
        nz = sum(1 for a in range(L.NA) if av[a])
        if s > out['score']:
            out['score'], out['nz'], out['sign'] = s, nz, (1 if d == delta else -1)
            out['czero'] = (av[c] == 0)
        if s > BASE:
            T.save(v, os.path.join(HERE, 'RJ_%d_c%d_u%d.json' % (s, c, u)))
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(CAND)
sweep(tag, CAND, evaluate, start, min(end, len(CAND)),
      keyfn=lambda s: '%d_%d' % s, budget=480)
rs = load(tag)
if rs:
    rs.sort(key=lambda r: -r['score'])
    print(f'\nbest residue jumps from {BASE}:')
    for r in rs[:12]:
        print('   a%-6d x%-6d -> %d  (nonzero atoms %s, target zeroed %s)'
              % (r['c'], r['u'], r['score'], r.get('nz'), r.get('czero')))
