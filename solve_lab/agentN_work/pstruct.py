"""Does the p-obstruction survive the polynomial framing, and what marks the nine carrier atoms?

Part 1 — the obstruction.  `obstruct.py` found that in all 924 six-row subsets the integrality
obstruction's denominator is divisible by
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663 .
That was measured on the NARROW 7-knob affine model.  Here the same question is asked of the exact
polynomial system's rank-14 variety built from the complete 68-knob set: reduce the 12x14 system
mod p and see whether the 7 unzeroable constants are exactly the ones that fall outside the mod-p
column space.

Part 2 — the carriers.  For the nine atoms that survive every re-orientation, measure the two
properties the polynomial framing can see: the exact degree of the atom as a polynomial in the
region's knobs, and whether the atom admits a legal unit target (i.e. can ever be a definition).
"""
import os, sys, json, pickle, ast, re
from collections import defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, model, optN
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits
from polyexact import P
from polyfull import exact_polys

Pp = 115792089237316195423570985008687907853269984665640564039457584007908834671663
CARRIERS = [22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887]

D = pickle.load(open(os.path.join(HERE, 'runs', 'solve68.pkl'), 'rb'))
M, b, Rl = D['M'], D['b'], D['Rl']
n = len(M[0])

print('=== PART 1: the p-obstruction on the exact polynomial variety ===', flush=True)
print('p is prime-ish check: p mod small primes ...', flush=True)
print('system %d x %d ; constants nonzero in %d rows' % (len(M), n, sum(1 for x in b if x)),
      flush=True)


def rank_mod(rows, ncol, p):
    A = [[x % p for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(ncol)]
        r += 1
    return r


def rank_q(rows, ncol):
    A = [[Fraction(x) for x in r] for r in rows]
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(A)):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(ncol)]
        r += 1
    return r


rQ = rank_q(M, n)
rP = rank_mod(M, n, Pp)
print('rank of the 12x%d region matrix : over Q = %d ; mod p = %d  (drop %d)'
      % (n, rQ, rP, rQ - rP), flush=True)
aug = [M[i] + [b[i]] for i in range(len(M))]
print('rank of [M|b]                   : over Q = %d ; mod p = %d'
      % (rank_q(aug, n + 1), rank_mod(aug, n + 1, Pp)), flush=True)
print('b mod p is zero in rows: %s'
      % [Rl[i] for i in range(len(M)) if b[i] % Pp == 0], flush=True)
print('b mod p nonzero in rows: %s'
      % [Rl[i] for i in range(len(M)) if b[i] % Pp], flush=True)
mzero = [Rl[i] for i in range(len(M)) if all(x % Pp == 0 for x in M[i])]
print('rows whose ENTIRE knob response vanishes mod p: %d -> %s' % (len(mzero), mzero), flush=True)
print('   (such a row demands b_i = 0 mod p, which no knob setting can change)', flush=True)
stuck = [Rl[i] for i in range(len(M)) if all(x % Pp == 0 for x in M[i]) and b[i] % Pp]
print('   of those, rows with b_i != 0 mod p — UNZEROABLE for a mod-p reason: %d -> %s'
      % (len(stuck), stuck), flush=True)

print('\n=== PART 2: what marks the nine carrier atoms ===', flush=True)
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']


def unit_targets(s):
    """variables that could legally be made the definition target of this atom:
    a bare `x_v` term at top level with coefficient +-1."""
    t = ast.parse(s, mode='eval').body
    out = set()

    def walk(node, sgn):
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
            walk(node.left, sgn)
            walk(node.right, sgn if isinstance(node.op, ast.Add) else -sgn)
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            walk(node.operand, -sgn)
        elif isinstance(node, ast.Name):
            out.add(int(node.id[2:]))
    walk(t, 1)
    return out


st = make(WIT)
b0 = build(st)
Rgn = b0['R']
atoms_R = set()
for e in Rgn:
    for c, a in ev.eq_terms[e][2]:
        atoms_R.add(a)
cands = set()
for q in atoms_R:
    if q in fr.csup:
        cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
cands = sorted(y for y in cands if y in FREE)

# exact polynomial of each atom of the region, over the 68 knobs
P.NK = len(cands)
import frameB as FB
v = list(st.v)
ns = {'v': v, '__builtins__': {}}
aff, ck = set(), set()
for j, Y in enumerate(cands):
    v[Y] = P.var(j, st.fv.get(Y, 0))
    aff.update(fr.desc[Y])
    ck.update(fr.chk[Y])
for u in sorted(aff, key=lambda u: fr.pos[u]):
    v[u] = eval(FB.DEFEXPR[u], ns)
apoly = {}
for a in sorted(ck):
    x = eval(FB.ACODE[a], ns)
    apoly[a] = x

print('%-8s %-6s %-8s %-7s %-6s %-9s %s' %
      ('atom', 'inR', 'nzWit', 'deg', '#eqs', '#unitTgt', 'note'), flush=True)
rowsout = []
for a in sorted(set(CARRIERS) | set(list(atoms_R)[:0])):
    inR = a in atoms_R
    x = apoly.get(a)
    deg = x.deg() if isinstance(x, P) else (0 if a not in apoly else -1)
    nz = st.av.get(a, 0) != 0
    ut = unit_targets(atom_src[a])
    neq = len(atom_eqs[a])
    print('%-8d %-6s %-8s %-7s %-6d %-9d %s'
          % (a, inR, nz, deg if isinstance(x, P) else 'const', neq, len(ut),
             'SQUARE (no bare var)' if not ut else ''), flush=True)
    rowsout.append(dict(atom=a, in_region=bool(inR), nonzero_at_witness=bool(nz),
                        poly_degree=(deg if isinstance(x, P) else None),
                        equations=neq, unit_targets=len(ut)))

# comparison population: every check atom of the region
print('\ncomparison over ALL %d atoms occurring in the region:' % len(atoms_R), flush=True)
dd = defaultdict(int)
ut0 = 0
for a in atoms_R:
    x = apoly.get(a)
    dd[x.deg() if isinstance(x, P) else 'const'] += 1
    if not unit_targets(atom_src[a]):
        ut0 += 1
print('   degree distribution: %s' % dict(dd), flush=True)
print('   atoms with NO legal unit target: %d of %d' % (ut0, len(atoms_R)), flush=True)

json.dump(dict(rank_Q=rQ, rank_modp=rP, stuck=stuck, mzero=mzero, carriers=rowsout),
          open(os.path.join(HERE, 'runs', 'pstruct.json'), 'w'), indent=1)
print('\nwrote runs/pstruct.json', flush=True)
