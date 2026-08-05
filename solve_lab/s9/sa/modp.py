"""Mod-P screening of the exact linear model (P = secp256k1 prime, prime => field)."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

P = lib.P


def solvable(M, b, n):
    """Is M x = b solvable over F_P?  M rows are lists of ints."""
    rows = [[c % P for c in r] + [bb % P] for r, bb in zip(M, b)]
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                pr = i; break
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        inv = pow(rows[r][c], P - 2, P)
        rows[r] = [(x * inv) % P for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(a - f * bq) % P for a, bq in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows):
            break
    for row in rows:
        if all(row[c] == 0 for c in range(n)) and row[n] != 0:
            return False
    return True
