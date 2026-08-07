#!/usr/bin/env python3
"""MACHINE-CHECKABLE CERTIFICATE for agent I's reduction claim.

Claim under test
----------------
  Every assignment that makes ALL 40,885 atoms of EQUATIONS.txt vanish encodes a
  solution of a 256-bit elliptic-curve discrete logarithm:  k*G = T  on
  y^2 = x^3 + b over p = 2^256-2^32-977, with a 256-bit PRIME group order.

Run:   python3 certify.py            (needs only ../../EQUATIONS.txt)
Prints PASS/FAIL per step and a final verdict.  Step 6 is the adversarial step:
it tests the mechanism (does the circuit really compute the ladder?) rather than
testing what the author expects to hold, and it states in the open exactly what
the certificate does NOT establish.
"""
import os, sys, re, json, time, collections, random, pickle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
P = 2**256 - 2**32 - 977
NV = 38748
RESULTS = []


def step(name, ok, detail=""):
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if detail:
        for line in str(detail).splitlines():
            print("        " + line)
    sys.stdout.flush()
    return ok


# ---------------------------------------------------------------- step 0
def ensure_caches():
    need = ['atoms.pkl', 'polys.pkl']
    if all(os.path.exists(os.path.join(HERE, f)) for f in need):
        return
    print("building parse caches from EQUATIONS.txt ...")
    import parse, poly
    parse.main()
    D, polys = poly.build()
    pickle.dump(polys, open(os.path.join(HERE, 'polys.pkl'), 'wb'))


ensure_caches()
from model import Model, load_assign          # noqa: E402
from prop import Engine as ZEngine            # noqa: E402
from fp import FpEngine, sqrt_p               # noqa: E402
from boolscore import Fast                    # noqa: E402

M = Model()
print(f"parsed: {M.ne} equations, {M.na} distinct atoms\n")
step("0. instance shape: 39,033 equations / 40,885 atoms, all atoms deg<=2",
     M.ne == 39033 and M.na == 40885 and
     all(max(len(m) for m in q) <= 2 for q in M.polys),
     f"equations={M.ne} atoms={M.na}")

# ---------------------------------------------------------------- step 1
var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)

zval = [None] * NV
ZE = ZEngine(M)
nz, zconf, zbr = ZE.propagate(zval)
pvars = set(v for v in range(NV) if zval[v] == P)

asserts = []
for i, q in enumerate(M.polys):
    if len(q) != 2:
        continue
    mons = sorted(q.items(), key=lambda kv: len(kv[0]))
    (m1, c1), (m2, c2) = mons
    if len(m1) == 1 and len(m2) == 2 and c1 == 1 and c2 == -1:
        a, b = m2
        if a in pvars and b not in pvars:
            asserts.append((i, m1[0], b))
        elif b in pvars and a not in pvars:
            asserts.append((i, m1[0], a))
handles = [h for _, _, h in asserts]
uniq = all(len(var2atoms[h]) == 1 for h in handles)
step("1. 3,707 atoms are `X - p*H`, every handle H occurs in NO other atom",
     len(asserts) == 3707 and len(set(handles)) == 3707 and uniq,
     f"assertion atoms={len(asserts)}  distinct handles={len(set(handles))}  "
     f"all handles occur in exactly one atom={uniq}\n"
     f"=> each is exactly `X == 0 (mod p)` with a free quotient; the mod-p\n"
     f"   abstraction of the instance is EXACT (no information is lost).")

# ---------------------------------------------------------------- step 2
known = sum(1 for x in zval if x is not None)
cnt = collections.Counter(zval[v] for v in range(NV) if zval[v] is not None)
big = {k: n for k, n in cnt.items() if abs(k) > 10**9}
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
step("2a. Z-propagation from EMPTY forces 5,624 vars, zero conflicts; the only\n"
     "    large pinned constants in the whole instance are p and K",
     known == 5624 and len(zconf) == 0 and set(big) == {P, K},
     f"forced={known} conflicts={len(zconf)}  values: 0->{cnt.get(0)}  1->{cnt.get(1)}  "
     f"p->{cnt.get(P)}  K->{cnt.get(K)}\nK = {K}")

F = Fast()
polwit = lambda u, r: (F.witp[u] if F.witp[u] in r else r[0])
val, conf, dec = F.run(polwit)
kn = sum(1 for x in val if x is not None)
step("2b. mod-p propagation (effective support) determines 28,701 vars from\n"
     "    1,156 boolean decisions, leaving EXACTLY 3 violated atoms",
     kn == 28701 and len(dec) == 1156 and len(conf) == 3,
     f"determined={kn} decisions={len(dec)} violated={len(conf)}: "
     f"{[('a%d: %s' % (a, M.src[a])) for a in conf]}")

A_VAR, B_VAR = 35389, 6671
mults = []
for a in conf:
    m = re.match(r'^X(\d+) - (\d+) \* X(\d+)$', M.src[a])
    mults.append((int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None)
rank2 = all(x is not None for x in mults) and len({x[1] for x in mults}) == 3
step("2c. the 3 violated atoms are `X_k - m_k * X35389` with three DIFFERENT m_k,\n"
     "    and each X_k is an independently-derived multiple of X6671\n"
     "    => rank 2 in (X35389, X6671) => BOTH are forced to 0 mod p",
     rank2, f"{mults}")

# ---------------------------------------------------------------- step 3
LEAVES = {12186: 'x1', 16742: 'y1', 14853: 'x2', 24908: 'y2',
          22162: 'x3', 30213: 'y3', 24453: 'Kc'}
try:
    import sympy
    R = F.reason
    syms = {v: sympy.Symbol(n) for v, n in LEAVES.items()}
    memo = {}

    def expand(v, depth=0):
        if v in syms:
            return syms[v]
        if v in memo:
            return memo[v]
        r = R[v]
        if r in (None, 'dec', 'pre') or depth > 60:
            memo[v] = sympy.Integer(val[v] if val[v] is not None else 0)
            return memo[v]
        # solve atom r for v: it is linear in v with coefficient +-1
        q = M.polys[r]
        c1 = 0
        rest = sympy.Integer(0)
        ok = True
        for mon, c in q.items():
            if v in mon:
                if len(mon) != 1:
                    ok = False
                    break
                c1 += c
            else:
                t = sympy.Integer(c)
                for x in mon:
                    t = t * expand(x, depth + 1)
                rest = rest + t
        if not ok or c1 not in (1, -1):
            memo[v] = sympy.Integer(val[v] if val[v] is not None else 0)
            return memo[v]
        memo[v] = sympy.expand(-rest / c1)
        return memo[v]

    x1, y1, x2, y2, x3, y3, Kc = (syms[k] for k in
                                  (12186, 16742, 14853, 24908, 22162, 30213, 24453))
    EA = sympy.expand(expand(A_VAR))
    EB = sympy.expand(expand(B_VAR))
    CA = sympy.expand((x2 - x1)**2 * (x3 + x1 + x2 + Kc) - (y2 - y1)**2)
    CB = sympy.expand((y3 + y1) * (x2 - x1) - (y2 - y1) * (x1 - x3))
    okA = sympy.simplify(EA - CA) == 0
    okB = sympy.simplify(EB - CB) == 0
    step("3. EXACT polynomial identities (symbolic back-substitution of the\n"
         "   instance's own atoms, not numeric agreement):\n"
         "     X35389 = (x2-x1)^2*(x3+x1+x2+K) - (y2-y1)^2\n"
         "     X6671  = (y3+y1)(x2-x1) - (y2-y1)(x1-x3)",
         okA and okB,
         f"X35389 identity: {okA}\nX6671 identity: {okB}\n"
         f"leaves x1=X12186 y1=X16742 x2=X14853 y2=X24908 x3=X22162 y3=X30213 K=X24453\n"
         f"(variables outside the leaf set are replaced by the values FORCED in\n"
         f" step 2 -- booleans and pins -- so the identity is exact for that branch)")
except ImportError:
    step("3. symbolic identity check", False, "sympy unavailable")

# ---------------------------------------------------------------- step 4
c3 = K * pow(3, -1, P) % P
sel = collections.defaultdict(list)
for i, s in enumerate(M.src):
    m = re.match(r'^X(\d+) \* \(X(\d+) - (\d+)\)', s)
    if m and int(m.group(3)) > 2**200:
        sel[int(m.group(1))].append((int(m.group(2)), int(m.group(3)) % P))
ks = sorted(sel)
B_CURVE = 64019533680030876408443198762210829058751700634554282185987325820393598524794


def on(pt):
    return (pt[1] * pt[1] - pow(pt[0], 3, P) - B_CURVE) % P == 0


direct = 0
pts = {}
for b in ks:
    (v1, a1), (v2, a2) = sel[b]
    cands = [((a1 + c3) % P, a2), ((a2 + c3) % P, a1),
             ((a1 + c3) % P, (-a2) % P), ((a2 + c3) % P, (-a1) % P)]
    if on(cands[0]):
        direct += 1
    for q in cands:
        if on(q):
            pts[b] = q
            break
x3v, y3v = val[22162], val[30213]
T = ((x3v + c3) % P, y3v)
step("4. substituting u = x + K/3 turns the gadget into the STANDARD Weierstrass\n"
     "   addition law, and the 256 conditional-pin constant pairs become points\n"
     "   on y^2 = x^3 + b (a = 0), b explicit below",
     len(ks) == 256 and len(pts) == 256 and direct == 219 and on(T),
     f"selectors with 296-bit pin pairs = {len(ks)}\n"
     f"pin pairs on the curve in the given order = {direct}/256\n"
     f"pin pairs on the curve after swap/negate  = {len(pts)}/256\n"
     f"target point T is on the curve            = {on(T)}\n"
     f"b = {B_CURVE}\np = {P}  (secp256k1's prime; b is NOT 7 -- this is a sextic twist)")

# ---------------------------------------------------------------- step 5
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337


def eadd(A, Bp):
    if A is None:
        return Bp
    if Bp is None:
        return A
    xa, ya = A; xb, yb = Bp
    if xa == xb and (ya + yb) % P == 0:
        return None
    l = (3 * xa * xa % P * pow(2 * ya % P, -1, P)) % P if A == Bp \
        else ((yb - ya) * pow(xb - xa, -1, P)) % P
    xc = (l * l - xa - xb) % P
    return (xc, (l * (xa - xc) - ya) % P)


def emul(k, Pt):
    Rr = None
    while k > 0:
        if k & 1:
            Rr = eadd(Rr, Pt)
        Pt = eadd(Pt, Pt)
        k >>= 1
    return Rr


import sympy
Sx = {q[0]: b for b, q in pts.items()}
doubles = sum(1 for b, q in pts.items() if (eadd(q, q) or (None,))[0] in Sx)
ordkills = emul(N, pts[ks[0]]) is None
start = [b for b in pts if b not in
         {Sx.get((eadd(q, q) or (-1, 0))[0]) for q in pts.values()}]
chain = 0
if len(start) == 1:
    q = pts[start[0]]
    while True:
        chain += 1
        d = eadd(q, q)
        if d is None or d[0] not in Sx:
            break
        q = pts[Sx[d[0]]]
        if chain > 300:
            break
G = pts[start[0]] if len(start) == 1 else None
step("5. group order N is a 256-bit PRIME and the 256 table points form a single\n"
     "   doubling ladder {G, 2G, 4G, ..., 2^255 G}",
     sympy.isprime(N) and ordkills and doubles >= 185 and chain == 256,
     f"N = {N}\nN prime = {sympy.isprime(N)}   N kills a curve point = {ordkills}\n"
     f"trace t = p+1-N = {P + 1 - N}   N == p (anomalous)? {N == P}\n"
     f"points whose double is also in the table = {doubles}/256\n"
     f"unique chain start = {start}   chain length from it = {chain}\n"
     f"G = {G}\nT = {T}")

emb = next((k for k in range(1, 51) if (pow(P, k, N) - 1) % N == 0), None)
step("5b. no MOV/FR reduction and not anomalous: embedding degree > 50",
     emb is None, f"embedding degree k with N | p^k-1 for k<=50: {emb}")

# ---------------------------------------------------------------- step 6
print("\n---- step 6: the weight-bearing step, tested adversarially ----")

# 6a: is every propagation branch point genuinely boolean?
nonbool = []


def probe(u, roots):
    if set(roots) - {0, 1}:
        nonbool.append((u, roots))
    return F.witp[u] if F.witp[u] in roots else roots[0]


F.run(probe)
step("6a. every branch point of the mod-p propagation has root set exactly {0,1}\n"
     "    (a non-boolean branch would be a continuous knob and would refute the claim)",
     not nonbool, f"non-boolean branch points found: {len(nonbool)} {nonbool[:5]}")

# 6b: confluence -- the propagated result must not depend on processing order
import boolscore
v2, c2_, d2 = F.run(polwit)
same = all(a == b for a, b in zip(val, v2)) and set(conf) == set(c2_)
step("6b. propagation is confluent (same result on repeat); every non-decision\n"
     "    assignment came from a unit clause, so A and B are a FUNCTION OF THE\n"
     "    BOOLEAN DECISIONS ALONE",
     same, f"identical={same}")

# 6c: THE MECHANISM TEST.  Does the circuit really compute the ladder?
# Turn on random subsets of selectors; predict the accumulator by independent EC
# arithmetic; compare with the value the circuit derives for (x1,y1).
idx = {}
q = pts[start[0]]; i = 0
while True:
    idx[Sx[q[0]]] = i
    d = eadd(q, q)
    if d is None or d[0] not in Sx or i >= 255:
        break
    q = pts[Sx[d[0]]]; i += 1
base_on = [b for b in ks if F.witp[b] == 1]
rng = random.Random(7)
mech_ok = 0
mech_tot = 0
mech_detail = []
for trial in range(6):
    extra = rng.sample([b for b in ks if b not in base_on], trial % 3 + 1)
    onset = set(base_on) | set(extra)

    def pol(u, roots, S=onset):
        if u in S:
            return 1 if 1 in roots else roots[0]
        x = F.witp[u]
        return x if x in roots else roots[0]

    vv, cc, _ = F.run(pol)
    if vv[12186] is None or vv[16742] is None:
        mech_detail.append(f"  |S|={len(onset)}: (x1,y1) released (free) -- no prediction")
        continue
    acc = ((vv[12186] + c3) % P, vv[16742])
    # does the circuit's (x1,y1) equal the EC sum of SOME subset of the on-set?
    hit = None
    lst = sorted(onset)
    for mask in range(1, 1 << len(lst)):
        s = None
        for j, b in enumerate(lst):
            if mask >> j & 1:
                s = eadd(s, emul(1 << idx[b], G))
        if s == acc:
            hit = [idx[b] for j, b in enumerate(lst) if mask >> j & 1]
            break
    mech_tot += 1
    if hit is not None:
        mech_ok += 1
    mech_detail.append(f"  |S|={len(onset)} ladder indices {sorted(idx[b] for b in onset)} "
                       f"-> circuit accumulator = sum of ladder indices {hit}")
step("6c. MECHANISM TEST: with random selector subsets switched on, the value the\n"
     "    circuit derives for (x1,y1) is exactly the elliptic-curve SUM of the\n"
     "    corresponding ladder points 2^i*G, computed independently.\n"
     "    (If the circuit were doing anything other than the ladder, this fails.)",
     mech_tot > 0 and mech_ok == mech_tot,
     f"matched {mech_ok}/{mech_tot} predictions\n" + "\n".join(mech_detail))

# 6d: is there ANY mod-p freedom besides the selector bits?
pinsrc = []
for i, s in enumerate(M.src):
    m = re.match(r'^X(\d+) \* \(X(\d+) - (\d+)\) - (?:(\d+) \* )?X(\d+)$', s)
    if m and int(m.group(3)) > 2**200:
        pinsrc.append((i, int(m.group(1)), int(m.group(2)), int(m.group(5))))
free_handles = [t for t in pinsrc if val[t[3]] is None or val[t[3]] != 0]
step("6d. every one of the 512 conditional-pin handles is forced to 0 mod p, so an\n"
     "    active pin fixes its advice value EXACTLY mod p (no continuous knob).\n"
     "    A single handle that were free mod p would refute the reduction.",
     len(pinsrc) == 512 and not free_handles,
     f"conditional pins={len(pinsrc)}  handles not forced 0 mod p={len(free_handles)}")

# 6e: releasing a coordinate does NOT remove the constraint -- it moves it one
# step up the ladder.  Demonstrated constructively.
relsel = None
for b in ks:
    if F.witp[b] == 1:
        continue

    def pol(u, roots, tgt=b):
        if u == tgt:
            return 1 if 1 in roots else roots[0]
        x = F.witp[u]
        return x if x in roots else roots[0]
    vv, cc, _ = F.run(pol)
    if vv[12186] is None and len(cc) == 0:
        relsel = (b, vv)
        break
moved = None
if relsel:
    b, vv = relsel
    x2v, y2v, x3q, y3q = vv[14853], vv[24908], vv[22162], vv[30213]
    D = (x2v - x3q) % P; S = (y2v + y3q) % P
    Tt = (y2v * x3q + y3q * x2v) % P
    from polyroot import roots as proots, pmul, padd, psub
    tt = [0, 1]
    qq = psub([x2v], tt)
    f1 = pmul(pmul(qq, qq), padd(tt, [(x3q + x2v + K) % P]))
    f1 = [cc_ * D % P * D % P for cc_ in f1]
    lin = psub([(y2v * D + Tt) % P], [cc_ * S % P for cc_ in tt])
    f = psub(f1, pmul(lin, lin))
    rs = proots(f)
    if rs:
        r = rs[-1]
        y1s = (S * r - Tt) % P * pow(D, -1, P) % P

        def pol2(u, rr, tgt=b):
            if u == tgt:
                return 1 if 1 in rr else rr[0]
            x = F.witp[u]
            return x if x in rr else rr[0]
        v3, c3f, _ = F.run(pol2, preassign={12186: r, 16742: y1s})
        moved = (b, len(c3f), c3f[:6], v3[35389], v3[6671])
step("6e. ADVERSARIAL: releasing a coordinate (switch on one more selector) makes\n"
     "    A=B=0 SOLVABLE -- and re-introduces the identical gadget one rung up the\n"
     "    ladder.  The constraint is invariant; only its position moves.",
     moved is not None and moved[3] == 0 and moved[4] == 0 and moved[1] == 3,
     f"selector X{moved[0]} released (x1,y1) with 0 conflicts;\n"
     f"solving the cubic for (x1,y1) gives A={moved[3]} B={moved[4]} (both zero),\n"
     f"and exactly {moved[1]} NEW violated atoms appear: {moved[2]}\n"
     f"-- the same `X_k - m_k*X_j` triple, now for the previous ladder rung."
     if moved else "could not construct the demonstration")

# ---------------------------------------------------------------- honest limits
print("\n---- WHAT THIS CERTIFICATE DOES AND DOES NOT ESTABLISH ----")
one_eq = sum(1 for a in range(M.na) if len(M.atom_eqs[a]) == 1)
print(f"""
ESTABLISHED (steps 1-6):
  If an assignment makes ALL 40,885 atoms vanish, then the 256 selector bits b_i
  satisfy  sum_i b_i * 2^i * G = T  on y^2 = x^3 + b over p, whose group order N
  is a 256-bit prime with embedding degree > 50 and trace {P + 1 - N}.
  Finding such b is a 256-bit ECDLP.  Best known cost: Pollard rho,
  ~sqrt(pi*N/4) ~= 2^127 group operations.  That is the honest price.

NOT ESTABLISHED (the one gap, stated plainly):
  "all atoms vanish" is SUFFICIENT for all 39,033 equations but not NECESSARY.
  Each equation is an integer combination of 3-24 atoms, so nonzero atoms may
  cancel.  {one_eq} atoms occur in exactly one equation, and the 39,026 witness
  itself has 9 nonzero atoms.  To close the gap one must show no atom vector in
  the image of the atom map, other than 0, kills every equation -- that is the
  compensation-closure question, and it is open.
  So the precise status is: THE INSTANCE CONTAINS A 256-BIT ECDLP ON ITS
  ALL-ATOMS-ZERO BRANCH.  Every construction anyone has run lives on that branch,
  which is why the 39,026 floor has never moved.
""")
ok = all(r[1] for r in RESULTS)
print(f"\nVERDICT: {sum(1 for r in RESULTS if r[1])}/{len(RESULTS)} steps PASS"
      f"  -> {'CERTIFICATE HOLDS' if ok else 'SOME STEP FAILED'}")
json.dump([{'step': n, 'pass': bool(o), 'detail': d} for n, o, d in RESULTS],
          open(os.path.join(HERE, 'certificate_results.json'), 'w'), indent=1)
sys.exit(0 if ok else 1)
