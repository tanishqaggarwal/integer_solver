"""U34: beam search over the DRV knobs, targeting the five discount equations.

Greedy on the failing count stalls at 11 because zeroing one of the five usually needs a
plateau move first.  Objective is lexicographic: (# of the five still non-zero, failing count).
Moves are exact integer roots of "equation e, as a polynomial in knob d", degree-verified.
usage: python3 u34_beam.py [beta a b]        (default: the ROOT with a generic honest leaf)
"""
import sys, math, time, json, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work/mirror')
import u20_sweep as S
import umodel as U, uscore as SC, checker
import harness as H
from u33_drvsolve import eqsum, solve_knob, DISC, DRV

ENG = SC.ENG
M = pickle.load(open('u_drvmap.pkl', 'rb'))['touch']       # knob -> equations it appears in
PAIRS = [(d, e) for d in DRV for e in DISC if e in M.get(d, ())]
BEAM, DEPTH = 4, 7


def state_of(seed):
    v = ENG.forward(seed)
    nz = tuple(e for e in DISC if eqsum(v, e) != 0)
    return (len(nz), len(checker.evaluate_all(SC.CODES, v))), nz, v


def beam(seed0, tag, tlimit=600):
    t0 = time.time()
    key0, nz0, _ = state_of(seed0)
    frontier = [(key0, seed0)]
    best = (key0, seed0)
    seen = {tuple(sorted(seed0.items()))}
    for lvl in range(DEPTH):
        cand = []
        for key, seed in frontier:
            for d, e in PAIRS:
                if time.time() - t0 > tlimit:
                    break
                for root in solve_knob(seed, e, d):
                    s2 = dict(seed)
                    if root == 0:
                        s2.pop(d, None)
                    else:
                        s2[d] = root
                    sig = tuple(sorted(s2.items()))
                    if sig in seen:
                        continue
                    seen.add(sig)
                    k2, nz2, _ = state_of(s2)
                    cand.append((k2, s2))
                    if k2 < best[0]:
                        best = (k2, s2)
                        print('   [%s] lvl%d x_%-6d:=root(eq%-5d) -> nz=%d failing=%d (%.0fs)'
                              % (tag, lvl, d, e, k2[0], k2[1], time.time() - t0))
                        sys.stdout.flush()
        if not cand or time.time() - t0 > tlimit:
            break
        cand.sort(key=lambda x: x[0])
        frontier = cand[:BEAM]
    return best


if __name__ == '__main__':
    if len(sys.argv) > 3:
        beta, a, b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    else:
        beta, a, b = U.ROOT, 47, 2081
    v, isl, valn = S.build(beta, a, b, a)
    for tag, extra in (('delivDRV', S.DRVSEED), ('noDRV', None)):
        sd = SC.seed_of_build(v, extra)
        (nz, n), sol = beam(sd, '%d/%s' % (beta, tag))
        print(' beta=%-6d a=%-6d b=%-6d %-9s -> %d failing (%d of 5 discount eqs still nonzero)'
              % (beta, a, b, tag, n, nz))
        if n <= 11:
            vv = ENG.forward(sol)
            p = 'u_beam_%d_%d_%d.json' % (beta, a, n)
            json.dump({("x_%d" % i): vv[i] for i in range(38748) if vv[i] != 0}, open(p, 'w'))
            print(' dumped %s' % p)
