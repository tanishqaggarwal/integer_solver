"""Stage 2: from a full11 (1,1) pair solution whose residual is the a-fold trio
   {20649, 20652, 32148}, run the FULL non-boolean knob closure and try to repair it.
   The deliverable repairs exactly this trio, and pays 7 equations for it; anything
   cheaper beats 39,026."""
import sys, os, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, sparse, full11 as F
OD = '/home/user/integer_solver/solve_lab/agentO_work'
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663

a = int(sys.argv[1]); b = int(sys.argv[2])
maxr = int(sys.argv[3]) if len(sys.argv) > 3 else 6
maxv = int(sys.argv[4]) if len(sys.argv) > 4 else 8000

r = F.solve_pair(a, b, verbose=False)
n, ns, av = r
print(f'stage1 pair ({a},{b}): fails={n} score={39033-n} bad={av}', flush=True)
v0 = E.forward(ns); bad0 = E.badatoms(v0)
S, cols, nonlin, rounds = simO.closure(v0, bad0, {18956, a, b}, maxr, maxv, verbose=True)
print('knobs', len(S), flush=True)
for at in sorted(bad0):
    reach = {f: cols[f][at] for f in S if at in cols[f]}
    lin = {f: c for f, c in reach.items() if (f, at) not in nonlin}
    print(f'row {at}: "{H.atoms[at][:70]}" reach={len(reach)} lin={len(lin)} '
          f'rhs_bits={abs(bad0[at]).bit_length()}', flush=True)
    for f, c in sorted(lin.items())[:12]:
        print(f'    x_{f}: ({len(str(abs(c)))}d) p|c={c%P==0}', flush=True)

# lazy repair
import lazy
t0 = time.time()
res = lazy.run(ns, {18956, a, b}, maxr=maxr, maxv=maxv, iters=30)
print(f'stage2 BEST fails={res[0]} score={39033-res[0]} bad={res[2]} ({time.time()-t0:.0f}s)', flush=True)
if res[0] < 7:
    json.dump({f"x_{i}": str(int(res[3][i])) for i in range(E.NV) if res[3][i] != 0},
              open(f'{OD}/stage2_{a}_{b}_{39033-res[0]}.json', 'w'))
    print('*** WROTE improvement', flush=True)
