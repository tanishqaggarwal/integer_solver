"""Measure the size and SHAPE of the irreducible residual problem.

Question: after all propagation, what is left?  Specifically
  (a) how many binary control bits actually move the two residues that must vanish,
  (b) what is the rank over GF(p) of the bit -> residue map (is there a subset-sum to solve?),
  (c) how many non-binary (256-bit integer) unknowns remain,
  (d) the algebraic degree of the residual conditions.
"""
import pickle, sys, collections
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
NV = 38748
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
freeinp = [x for x in range(NV) if x not in definer]
bfree = [b for b in freeinp if b in boolv]


def resid(v):
    """the two scalar conditions that must vanish for the whole system to hold"""
    D1 = (v[7068] - v[2099] - 7376877 * v[642]) % P
    D2 = (v[4432] - v[19964] - v[28730]) % P
    return D1, D2


if __name__ == '__main__':
    v0 = H.load_assignment('../best/new_instance_partial_39022.json')
    b0 = resid(v0)
    print(f'residual conditions at the partial:  D1%p = {b0[0]}\n'
          f'                                      D2%p = {b0[1]}')
    print(f'\nboolean free inputs: {len(bfree)}')
    movers = []
    inert = 0
    for b in bfree:
        v = list(v0); ch, _ = ripple(v, {b: 1 - v0[b]})
        d = resid(v)
        if d != b0:
            movers.append((b, ((d[0]-b0[0]) % P, (d[1]-b0[1]) % P)))
        if not ch: inert += 1
    print(f'bits that move (D1,D2) mod p: {len(movers)}')
    for b, d in movers: print(f'   x_{b}: dD1={d[0]}  dD2={d[1]}')
    # rank over GF(p) of the bit->residue map
    rows = [list(d) for _, d in movers]
    rank = 0
    M = [r[:] for r in rows]
    for c in range(2):
        piv = next((i for i in range(rank, len(M)) if M[i][c] % P), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][c], P-2, P)
        for i in range(len(M)):
            if i != rank and M[i][c] % P:
                f = M[i][c] * inv % P
                for j in range(2): M[i][j] = (M[i][j] - f*M[rank][j]) % P
        rank += 1
    print(f'\nRANK of the bit -> (D1,D2) map over GF(p): {rank} / 2')
    # non-binary free inputs that move the residues
    nb = []
    for f in freeinp:
        if f in boolv: continue
        v = list(v0); ripple(v, {f: v0[f] + 1})
        if resid(v) != b0: nb.append(f)
    print(f'non-boolean free inputs that move (D1,D2) mod p: {len(nb)} -> {nb}')
