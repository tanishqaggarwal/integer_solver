"""S11 step 69: WHY the mod-p model mispredicts, and how to make it exact.

656 of 1,376 large-move predictions were wrong, which is what makes every linear
veto in this lab untrustworthy.  The cause is not the mod-p reduction -- it is
INTEGRALITY.  Forward evaluation solves each gate atom for its output over Z via
T.solve_lin, which returns None when the solution is not an integer; the variable
then keeps its stale value and the gate silently breaks.  Mod p there is no such
thing as a non-integral solution, so the model and the machine part ways.

Fix: move along a SUBLATTICE.  If every free input moves by a multiple of N, then
every atom's increment is a multiple of N (sums and products of multiples of N are
multiples of N), so dividing by any pivot c with c | N keeps the quotient integral
-- as long as N carries enough copies of each pivot for the DAG's depth.  Take
    D = lcm of every pivot seen in one forward pass,   N = D**k
and the model should become exact.  Mod p, N is invertible, so the sublattice hits
every residue: nothing is lost, only fidelity is gained.

Usage: fidelity.py [K] [NPROBE]
"""
import os, sys, random, time, math, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
random.seed(7)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NP = int(sys.argv[2]) if len(sys.argv) > 2 else 3

v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = sorted(a for a in range(L.NA) if a not in L.atom_out)
BAD = [21617, 29539]
FORBID = {2081, 4287}
U = sorted((set(ad.grad(BAD[0], vm0)) | set(ad.grad(BAD[1], vm0))) - FORBID)
cols = {u: jac_column(u, v0, vm0, CHECKS) for u in U}

# ---- the pivot lattice ------------------------------------------------------
piv = set()
for t in ad.ORDER:
    a = L.definer[t]
    d = ad.dpart(a, t, v0)          # exact integer pivot at the base state
    if d:
        piv.add(abs(d))
piv.discard(1)
small = sorted(x for x in piv if x < 10 ** 6)
print(f'{len(piv)} distinct pivots; {len(small)} below 10^6; '
      f'largest {max(piv):.3g}' if piv else 'no pivots', flush=True)
D = 1
for x in small:
    D = D * x // math.gcd(D, x)
# the huge pivots are variable-valued (quadratic gates); fold in their prime-ish
# content cheaply by multiplying, not lcm-ing, the ones we can afford
print(f'lcm of the small pivots: {D.bit_length()} bits', flush=True)
N = pow(D, K)
print(f'N = D**{K}: {N.bit_length()} bits;  gcd(N,p) = {math.gcd(N, P)}', flush=True)


def measure(step_fn, label):
    ok = collections.Counter()
    bad = collections.Counter()
    gates = 0
    for u in U:
        col = cols[u]
        for _ in range(NP):
            d = step_fn()
            w = list(v0)
            w[u] = w[u] + d
            ad.fwd(w, rounds=6)
            aw = L.all_atom_values(w)
            gates += sum(1 for a in L.atom_out if aw[a])
            for c in col:
                if (av0[c] + col[c] * d) % P == aw[c] % P:
                    ok[c] += 1
                else:
                    bad[c] += 1
    tot = sum(ok.values()) + sum(bad.values())
    rows = set(ok) | set(bad)
    exact = sum(1 for c in rows if bad[c] == 0)
    print(f'{label:28s} predictions {sum(ok.values())}/{tot} correct '
          f'({100.0*sum(ok.values())/max(tot,1):.1f}%);  rows exact '
          f'{exact}/{len(rows)};  broken gates/probe '
          f'{gates/(len(U)*NP):.1f}', flush=True)
    return exact, len(rows)


t0 = time.time()
measure(lambda: random.randrange(1, P), 'free (any integer)')
measure(lambda: N * random.randrange(1, 1 << 40), f'sublattice N = D**{K}')
measure(lambda: D * random.randrange(1, 1 << 40), 'sublattice N = D')
print(f'({time.time()-t0:.0f}s)')
