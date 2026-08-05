"""Maximise how many of the 13 sacrificed equations vanish, over the exact integer knob lattice.

Knobs: z = (t, h, s, k, y5, y6, y7, y8)
    y1 = atom 22229 = D1 - 7376877*t        (t   = x_642)
    y2 = atom 35762 = t - p*h               (h   = x_17325)
    y3 = atom 22230 = s - p*k               (s   = x_28730, k = x_9413)
    y4 = atom 22231 = D2 - s
    y5 = atom 35758, y6 = atom 35759        free via (x_29854, x_1329, x_9118)  [CRT: gcd(5113045,p)=1]
    y7 = atom 35760, y8 = atom 35761        free via (x_31864, x_10903, x_8731)
    y9 = atom 37887 = y4**2
Rational-rank prefilter, then exact Z solvability via Smith Normal Form.
"""
import pickle, itertools, sys, os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
S9 = os.path.dirname(HERE)
sys.path.insert(0, S9); sys.path.insert(0, HERE); os.chdir(S9)
import harness as H
from snf import solve_int

P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); eq_terms = d['eq_terms']
v = H.load_assignment('../best/new_instance_partial_39022.json')
D1 = v[7068] - v[2099] - 7376877 * v[642]
D2 = v[4432] - v[19964]
ATOMS = [22229, 35762, 22230, 22231, 35758, 35759, 35760, 35761, 37887]
S13 = [2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125, 9123, 18673]
NZ = 8


def ycol(j):
    m = [0] * 8
    if j == 0: m[0] = -7376877; m[1] = 1
    if j == 1: m[1] = -P
    if j == 2: m[2] = 1; m[3] = -1
    if j == 3: m[2] = -P
    if j >= 4: m[j] = 1
    return m


YC = [ycol(j) for j in range(NZ)]
CONST = [D1, 0, 0, D2, 0, 0, 0, 0]

co_by_eq = {}
for i in S13:
    m, sq, tl = eq_terms[i]
    co = {a: 0 for a in ATOMS}
    for c, a in tl:
        if a in co: co[a] += c
    co_by_eq[i] = co


def rows_for(T):
    A = []; b = []
    for i in T:
        co = dict(co_by_eq[i])
        if co.get(37887, 0):
            A.append([1 if j == 2 else 0 for j in range(NZ)]); b.append(D2)
            co[37887] = 0
        row = [0] * NZ; rhs = 0
        for a, c in co.items():
            if not c: continue
            kk = ATOMS.index(a)
            if kk >= 8: continue
            for j in range(NZ): row[j] += c * YC[j][kk]
            rhs -= c * CONST[kk]
        if any(row): A.append(row); b.append(rhs)
    return A, b


def q_feasible(A, b):
    if not A: return True
    m = len(A); n = len(A[0])
    M = [[Fraction(x) for x in A[k]] + [Fraction(b[k])] for k in range(m)]
    r = 0
    for c in range(n):
        piv = None
        for k in range(r, m):
            if M[k][c] != 0: piv = k; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        for k in range(m):
            if k != r and M[k][c] != 0:
                f = M[k][c] / pv
                for j in range(c, n + 1): M[k][j] -= f * M[r][j]
        r += 1
        if r == m: break
    for k in range(r, m):
        if M[k][n] != 0: return False
    return True


if __name__ == '__main__':
    known = [2554, 6816, 8124, 8680]
    A, b = rows_for(known)
    print('validation: known 4-subset  Q-feasible=%s  Z-solvable=%s'
          % (q_feasible(A, b), solve_int(A, b) is not None), flush=True)
    best = None
    for size in range(len(S13), 0, -1):
        hits = []
        for T in itertools.combinations(S13, size):
            A, b = rows_for(list(T))
            if q_feasible(A, b): hits.append(T)
        print(f'size {size}: {len(hits)} Q-feasible subsets', flush=True)
        for T in hits:
            A, b = rows_for(list(T))
            z = solve_int(A, b) if A else [0] * NZ
            if z is not None:
                best = (T, z)
                print(f'\n*** MAX zeroable = {size}: {T}\n    z = {z}'
                      f'\n    failing = {len(S13) - size} -> score {39033 - (len(S13) - size)}', flush=True)
                break
        if best: break
    pickle.dump(best, open('kernel/opt7.pkl', 'wb'))
