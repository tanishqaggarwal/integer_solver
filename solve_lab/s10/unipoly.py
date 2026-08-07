"""S11 step 78: EXACT univariate analysis -- the first non-linearised model of a move.

Every gate output coefficient in this instance is +-1 (30,418 of +1 and 1,057 of -1),
so forward evaluation performs NO division at all: the map from the free inputs to
every atom value is an honest polynomial over Z, hence over F_p.  That licenses
something the whole lab has been unable to do -- model a move EXACTLY.

For one free input u, evaluate the instance at x_u + d for d = 0..K.  Each check c
then gives K+1 samples of the univariate polynomial f_c(d) = c(v + d*e_u) mod p, and
Newton interpolation recovers f_c exactly (verified at extra points).  With the real
polynomials in hand:

  * checks that currently HOLD and vary with u have f_c(0) = 0, so x | f_c.  The gcd
    G of all of them, divided by x, has as its roots exactly the nonzero jumps along
    u that break NOTHING -- nonlinear symmetries invisible to every Jacobian here.
  * checks that currently FAIL need f_c(d) = 0.  gcd(G, those) has as its roots the
    jumps that fix them while breaking nothing.  A single root is a strict win.

Root finding is Cantor-Zassenhaus over F_p.  Nothing is linearised and nothing is
predicted: the polynomials are exact and the roots are exact.

Usage: unipoly.py START END [state.json] [K]
"""
import os, sys, time, random
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import suppfree
from chunk import sweep, load
P = ad.P
random.seed(3)

# ---------- polynomials over F_p, coefficient lists low-order first -------------
def trim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def pmul(f, g):
    if not f or not g:
        return []
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                r[i + j] = (r[i + j] + a * b) % P
    return trim(r)


def psub(f, g):
    n = max(len(f), len(g))
    return trim([( (f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % P
                 for i in range(n)])


def pmod(f, g):
    f = f[:]
    dg = len(g) - 1
    inv = pow(g[-1], -1, P)
    while len(f) - 1 >= dg and f:
        k = len(f) - 1 - dg
        c = f[-1] * inv % P
        for i in range(dg + 1):
            f[k + i] = (f[k + i] - c * g[i]) % P
        trim(f)
    return f


def pgcd(f, g):
    f, g = trim(f[:]), trim(g[:])
    while g:
        f, g = g, pmod(f, g)
    if f:
        inv = pow(f[-1], -1, P)
        f = [x * inv % P for x in f]
    return f


def ppow(b, e, m):
    r, b = [1], pmod(b[:], m)
    while e:
        if e & 1:
            r = pmod(pmul(r, b), m)
        b = pmod(pmul(b, b), m)
        e >>= 1
    return r


def roots(f):
    """All roots in F_p of f (Cantor-Zassenhaus)."""
    f = trim(f[:])
    if len(f) <= 1:
        return []
    g = pgcd(psub(ppow([0, 1], P, f), [0, 1]), f)     # gcd(f, x^p - x)
    out = []
    st = [g]
    while st:
        h = st.pop()
        if len(h) <= 1:
            continue
        if len(h) == 2:
            out.append((-h[0]) * pow(h[1], -1, P) % P)
            continue
        for _ in range(40):
            a = random.randrange(P)
            t = pgcd(psub(ppow([a, 1], (P - 1) // 2, h), [1]), h)
            if 1 < len(t) < len(h):
                st.append(t)
                st.append(pmod(h, t) and h or h)      # placeholder, replaced below
                st.pop()
                q = h[:]
                # exact division h / t
                dq, inv = len(t) - 1, pow(t[-1], -1, P)
                res = [0] * (len(h) - len(t) + 1)
                while len(q) - 1 >= dq and q:
                    k = len(q) - 1 - dq
                    c = q[-1] * inv % P
                    res[k] = c
                    for i in range(dq + 1):
                        q[k + i] = (q[k + i] - c * t[i]) % P
                    trim(q)
                st.append(trim(res))
                break
        else:
            continue
    return sorted(set(out))


def interp(xs, ys):
    """Newton interpolation mod p; returns coefficients low-order first."""
    n = len(xs)
    dd = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) * pow(xs[i] - xs[i - j], -1, P) % P
    f = [0]
    for i in range(n - 1, -1, -1):
        f = psub(pmul(f, [(-xs[i]) % P, 1]), [(-dd[i]) % P])
    return trim(f)


# ---------- the instance -------------------------------------------------------
src = sys.argv[3] if len(sys.argv) > 3 else 'mod9118_0.json'
K = int(sys.argv[4]) if len(sys.argv) > 4 else 9
tag = 'uni_' + os.path.basename(src).replace('.json', '')
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
bav = L.all_atom_values(base)
BASE = L.NEQ - len(L.failing_eqs(bav))
BAD = [a for a in range(L.NA) if a not in L.atom_out and bav[a]]
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
print(f'{src}: score {BASE}; failing checks {BAD}; K = {K}', flush=True)
_, freelist, SVS = suppfree.build(base, modp=None)
US = set()
for c in BAD:
    m = suppfree.atom_supp(c, base, SVS, modp=None)
    US |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
US = sorted(US)
print(f'{len(US)} free inputs structurally reach a failing check', flush=True)


def evaluate(u):
    xs = list(range(K + 3))
    cols = []
    for d in xs:
        v = list(base)
        v[u] = v[u] + d
        ad.fwd(v, rounds=6)
        av = L.all_atom_values(v)
        cols.append([av[c] % P for c in CHECKS])
    vary = [i for i in range(len(CHECKS)) if len(set(c[i] for c in cols)) > 1]
    G, fails, deg_bad = None, [], 0
    for i in vary:
        ys = [cols[d][i] for d in xs[:K + 1]]
        f = interp(xs[:K + 1], ys)
        ok = all(sum(co * pow(x, e, P) for e, co in enumerate(f)) % P == cols[x][i]
                 for x in xs[K + 1:])
        if not ok:
            deg_bad += 1
            continue
        c = CHECKS[i]
        if bav[c] == 0:
            G = f if G is None else pgcd(G, f)
        else:
            fails.append(f)
    out = {'u': u, 'vary': len(vary), 'degfail': deg_bad,
           'Gdeg': (len(G) - 1) if G else None}
    if G is None:                       # nothing that holds varies with u: fully free
        out['safe'] = 'ALL'
        H = fails[0] if fails else None
        for f in fails[1:]:
            H = pgcd(H, f)
    else:
        H = G
        for f in fails:
            H = pgcd(H, f)
    if H:
        rs = [r for r in roots(H) if r % P]
        out['fixroots'] = len(rs)
        if rs:
            best = BASE
            for r in rs[:6]:
                for d in (r, r - P):
                    v = list(base)
                    v[u] = v[u] + d
                    ad.fwd(v, rounds=6)
                    av = L.all_atom_values(v)
                    s = L.NEQ - len(L.failing_eqs(av))
                    if s > best:
                        best = s
                        T.save(v, os.path.join(HERE, 'UNI_%d_u%d.json' % (s, u)))
            out['score'] = best
    if G:
        sr = [r for r in roots(G) if r % P]
        out['saferoots'] = len(sr)
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(US)
sweep(tag, US, evaluate, start, min(end, len(US)), keyfn=lambda u: str(u), budget=540)
rs = load(tag)
if rs:
    print('\nfree inputs with a nonzero root that breaks nothing:')
    for r in sorted(rs, key=lambda r: -(r.get('score') or 0))[:15]:
        print('   x%-6d varying checks %-4d  gcd deg %-4s  fix-roots %-3s  '
              'safe-roots %-3s  score %s'
              % (r['u'], r['vary'], r.get('Gdeg'), r.get('fixroots'),
                 r.get('saferoots'), r.get('score')))
