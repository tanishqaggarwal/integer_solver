"""S11 step 109: the deliverable's twelve equations are HOMOGENEOUS in the seven atoms.

At the 39,026 witness the twelve equations touched by the residual involve 24 atoms,
seventeen of which are zero, so each equation reduces to a linear form in the seven
residual values alpha -- with **zero constant term**:

    the 12 x 7 matrix has rank 7;  alpha = 0 satisfies all twelve
    the current alpha satisfies exactly five:  [2554, 6816, 8124, 9123, 9421]

And every alpha is EXACTLY realisable, because frame 2 detaches five of the outputs
and the rest are free:

    x642    = x17325*p + a6            x7068  = x2099 + 7376877*x642 + a0
    x28730  = x9413*p  + a1
    x29854  = 5113045*x7075*x9118 - a3     x1329  = (x29854 - a2)/p
    x31864  = a5 - x7075*x8731             x10903 = (x31864 - a4)/p

the last two needing only  5113045*x7075*x9118 = a2 + a3  and  x7075*x8731 = a5 - a4
(mod p), which fixes the RESIDUES of the two free inputs x9118 and x8731 and nothing
else.  So the seven atom values are free, and the only side effect of choosing them is
whatever those two residues do to the rest of the instance.

Since the system is homogeneous, zeroing k rows leaves a (7-k)-dimensional family:

    k = 7  ->  alpha = 0 only          build7's route, side-effect cost 29 -> 39,004
    k = 6  ->  a free scalar t         6 of 12 fail  ->  39,027 - cost   <-- BEATS 39,026
    k = 5  ->  2 parameters            7 of 12 fail  ->  39,026 - cost   <-- the deliverable

So the deliverable sits at the k = 5 optimum, and k = 6 is one better *if* some scalar
in the family has side-effect cost zero.  Enumerate the 6-subsets, take the kernel
vector, sweep the scalar, realise exactly, and score.

Usage: sixsub.py START END [NT]
"""
import os, sys, time, itertools, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd
from chunk import sweep, load
P = ad.P
random.seed(31)
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 6
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(av))
E = sorted(set(e for a in SEVEN for e in L.atom2eq[a]))
M = [[L.eq_atoms[e][2].get(a, 0) for a in SEVEN] for e in E]
print('witness score %d;  %d equations, homogeneous in the seven atoms'
      % (BASE, len(E)), flush=True)
cur = [av[a] for a in SEVEN]
SAT = [i for i in range(len(E)) if sum(M[i][j] * cur[j] for j in range(7)) == 0]
print('currently satisfied rows: %s (%d of %d)'
      % ([E[i] for i in SAT], len(SAT), len(E)), flush=True)


def kernel(rows):
    m = 7
    A = [r[:] for r in rows]
    n = len(A)
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, n) if A[i][j] % P), None)
        if k is None:
            continue
        A[r_], A[k] = A[k], A[r_]
        inv = pow(A[r_][j], -1, P)
        A[r_] = [x * inv % P for x in A[r_]]
        for i in range(n):
            if i != r_ and A[i][j] % P:
                f = A[i][j]
                A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
        piv.append(j)
        r_ += 1
    out = []
    for j0 in [j for j in range(m) if j not in set(piv)]:
        w = [0] * m
        w[j0] = 1
        for i, j in enumerate(piv):
            w[j] = (-A[i][j0]) % P
        out.append(w)
    return out, r_


def realise(alpha):
    """Write the seven atom values exactly.  Returns None if not realisable."""
    v = list(base)
    x7075 = v[7075] % P
    if x7075 % P == 0:
        return None
    # x9118 residue from  5113045*x7075*x9118 = a2 + a3
    c = 5113045 * x7075 % P
    r1 = (alpha[2] + alpha[3]) % P * pow(c, -1, P) % P
    # x8731 residue from  x7075*x8731 = a5 - a4
    r2 = (alpha[5] - alpha[4]) % P * pow(x7075, -1, P) % P
    v[9118] = (v[9118] // P) * P + r1
    v[8731] = (v[8731] // P) * P + r2
    fwd(v)
    x7075 = v[7075]
    n29854 = 5113045 * x7075 * v[9118] - alpha[3]
    if (n29854 - alpha[2]) % P:
        return None
    v[29854] = n29854
    v[1329] = (n29854 - alpha[2]) // P
    n31864 = alpha[5] - x7075 * v[8731]
    if (n31864 - alpha[4]) % P:
        return None
    v[31864] = n31864
    v[10903] = (n31864 - alpha[4]) // P
    v[28730] = v[9413] * P + alpha[1]
    v[642] = v[17325] * P + alpha[6]
    fwd(v)
    v[7068] = v[2099] + 7376877 * v[642] + alpha[0]
    fwd(v)
    return v


# sanity: reproduce the witness exactly
chk = realise(cur)
if chk is not None:
    print('sanity: re-realising the current alpha gives score %d (expect %d)'
          % (L.NEQ - len(L.failing_eqs(L.all_atom_values(chk))), BASE), flush=True)

SUBS = []
for S in itertools.combinations(range(len(E)), 6):
    ker, rk = kernel([M[i] for i in S])
    if rk == 6 and len(ker) == 1:
        SUBS.append((S, ker[0]))
print('%d six-subsets with a one-dimensional family' % len(SUBS), flush=True)
# put the ones containing all five currently-satisfied rows first
SS = set(SAT)
SUBS.sort(key=lambda sk: -len(SS & set(sk[0])))


def evaluate(spec):
    S, k = spec
    best, bt = -1, None
    for t in [1, 2, 3] + [random.randrange(1, P) for _ in range(NT)]:
        alpha = [t * x % P for x in k]
        v = realise(alpha)
        if v is None:
            continue
        aw = L.all_atom_values(v)
        s = L.NEQ - len(L.failing_eqs(aw))
        if s > best:
            best, bt = s, t
        if s > BASE:
            T.save(v, os.path.join(HERE, 'SIX_%d_%s.json'
                                   % (s, '_'.join(map(str, S)))))
    return {'S': list(S), 'overlap': len(SS & set(S)), 'score': best, 't': bt}


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(SUBS)
sweep('sixsub', SUBS, evaluate, start, min(end, len(SUBS)),
      keyfn=lambda sk: ','.join(map(str, sk[0])), budget=540)
rs = load('sixsub')
if rs:
    rs.sort(key=lambda r: -r['score'])
    print('\nbest six-subsets (witness is %d):' % BASE)
    for r in rs[:15]:
        print('   rows %-28s overlap %d  score %d'
              % ([E[i] for i in r['S']], r['overlap'], r['score']))
