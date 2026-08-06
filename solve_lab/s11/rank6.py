import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, sys6
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
BITSET = set(int(b) for b in LD)
S = json.load(open(os.path.join(HERE, 'data', 'sys6.json')))
NAMES = sys6.NAMES
pool = sorted({u for nm in NAMES for u in S[nm] if u not in BITSET})
print(f"non-bit control pool: {len(pool)}")
rnd = random.Random(3)
th = {u: sys6.BASE[u] + rnd.randrange(1, 1 << 40) for u in pool}
r = sys6.six(sys6.ev(th))
J = []
for i in range(6):
    J.append([0] * len(pool))
for j, c in enumerate(pool):
    t2 = dict(th)
    t2[c] = th[c] + 1
    r1 = sys6.six(sys6.ev(t2))
    for i in range(6):
        J[i][j] = (r1[i] - r[i]) % P


def rank_mod(Mx, p):
    A = [row[:] for row in Mx]
    m, n = len(A), len(A[0])
    rk = 0
    for c in range(n):
        pr = None
        for i in range(rk, m):
            if A[i][c] % p:
                pr = i
                break
        if pr is None:
            continue
        A[rk], A[pr] = A[pr], A[rk]
        inv = pow(A[rk][c], -1, p)
        A[rk] = [x * inv % p for x in A[rk]]
        for i in range(m):
            if i != rk and A[i][c] % p:
                f = A[i][c]
                A[i] = [(A[i][k] - f * A[rk][k]) % p for k in range(n)]
        rk += 1
        if rk == m:
            break
    return rk


print("rank(J)      =", rank_mod(J, P))
print("rank([J | r]) =", rank_mod([J[i] + [r[i]] for i in range(6)], P))
for i in range(6):
    nz = [pool[j] for j in range(len(pool)) if J[i][j]]
    print(f"  {NAMES[i]:8s}: {len(nz)} controls")
# which pairs of rows are proportional?
import itertools
for a, b in itertools.combinations(range(6), 2):
    ratio = None
    prop = True
    for j in range(len(pool)):
        x, y = J[a][j] % P, J[b][j] % P
        if x == 0 and y == 0:
            continue
        if x == 0 or y == 0:
            prop = False
            break
        rr = y * pow(x, -1, P) % P
        if ratio is None:
            ratio = rr
        elif rr != ratio:
            prop = False
            break
    if prop and ratio is not None:
        print(f"  ROWS {NAMES[a]} and {NAMES[b]} are PROPORTIONAL (factor {str(ratio)[:20]}...)")
