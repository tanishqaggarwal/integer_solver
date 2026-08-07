"""S11 step 77: the two-phase move -- residue jump (mod p) then integer lift (Z).

§124 splits the instance in two.  Mod p the handles are invisible (their coefficient
is a multiple of p), so the mod-p phase is the whole difficulty; over Z the handles
enter linearly with coefficient d*p, so once a check is ≡ 0 (mod p) zeroing it
exactly is arithmetic.  This is the first solver in the lab that does both:

  phase 1  residue jump:  x_u <- x_u + delta,  delta = -c/(dc/dx_u) mod p
           -- an exact jump of a full residue class, not a tangent step
  phase 2  integer lift:  for every check that is now ≡ 0 (mod p) but nonzero,
           find a free input whose EXACT INTEGER coefficient divides it and set it

Phase 2 needs intad.jacZ: ad.dpart reduces mod p and so reports 0 for every p-wire
multiplication, which is precisely where the handles live.  Every previous repair in
this lab was blind to them.

Usage: rj2.py START END [state.json]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
from intad import jacZ
import suppfree
from chunk import sweep, load
P = ad.P
src = sys.argv[3] if len(sys.argv) > 3 else 'mod9118_0.json'
tag = 'rj2_' + os.path.basename(src).replace('.json', '')

base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bm = [x % P for x in base]
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
BAD = [a for a in range(L.NA) if a not in L.atom_out and bav[a]]
print(f'{src}: score {BASE}; failing checks {BAD}', flush=True)

_, freelist, SVS = suppfree.build(base, modp=None)      # structural, value-free


def struct_supp(c, v):
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    return [freelist[i] for i in range(len(freelist)) if (m >> i) & 1]


def lift(v, av, budget=90):
    """Zero every check that is ≡0 mod p, using exact integer coefficients."""
    t0 = time.time()
    best = L.NEQ - len(L.failing_eqs(av))
    moved = True
    while moved and time.time() - t0 < budget:
        moved = False
        cand = [a for a in range(L.NA) if a not in L.atom_out and av[a]
                and av[a] % P == 0]
        for c in cand:
            for u in struct_supp(c, v):
                g = jacZ(u, v, [c]).get(c, 0)
                if not g or av[c] % g:
                    continue
                w = list(v)
                w[u] = w[u] - av[c] // g
                ad.fwd(w, rounds=6)
                aw = L.all_atom_values(w)
                s = L.NEQ - len(L.failing_eqs(aw))
                if aw[c] == 0 and s >= best:
                    v, av, best, moved = w, aw, s, True
                    break
            if moved:
                break
    return v, av, best


# candidates: every (failing check, free input that moves it mod p)
CAND = []
for c in BAD:
    for u in struct_supp(c, base):
        if jac_column(u, base, bm, [c]).get(c, 0):
            CAND.append((c, u))
print(f'{len(CAND)} residue jumps available', flush=True)


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
        pre = L.NEQ - len(L.failing_eqs(av))
        v, av, s = lift(v, av)
        if s > out['score']:
            out.update(score=s, pre=pre, nz=sum(1 for a in range(L.NA) if av[a]),
                       sign=1 if d == delta else -1, czero=(av[c] == 0))
        if s > BASE:
            T.save(v, os.path.join(HERE, 'RJ2_%d_c%d_u%d.json' % (s, c, u)))
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(CAND)
sweep(tag, CAND, evaluate, start, min(end, len(CAND)),
      keyfn=lambda s: '%d_%d' % s, budget=540)
rs = load(tag)
if rs:
    rs.sort(key=lambda r: -r['score'])
    print(f'\nbest two-phase moves from {BASE}:')
    for r in rs[:12]:
        print('   a%-6d x%-6d  jump %d -> lift %d  (nonzero atoms %s, target zeroed %s)'
              % (r['c'], r['u'], r.get('pre', -1), r['score'], r.get('nz'),
                 r.get('czero')))
