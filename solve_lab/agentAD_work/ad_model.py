"""Agent AD -- the small analogue of the EQUATIONS.txt merge-tree system.

Built from the DECODED MATHEMATICS ONLY (fleet rule 2): the law shapes recorded in
agentW_work/RESUME_W.md and the campaign brief.  Nothing here looks at the layout,
ordering or coefficient templates of EQUATIONS.txt.

------------------------------------------------------------------------------
THE LAW (real instance, per merge block; W's five-atom integer form)

    congruence k=1,2,3 :   a_k  * L     * ( c_k1*N1 + c_k2*N2 )  =  c_k * P * u_k
    off-pin    j=5,6   :   a'_j * (1-L) *   i_j                  =  c'_j* P * u'_j

    A  = i1 - i2      B  = i4 - i3       E = i1 + i2 + i5 + Q
    N1 = E*A^2 - B^2  N2 = A*(i3+i6) - B*(i2-i5)

Identifying the slots with curve coordinates (this is forced, and is checked in
ad_validate.py by evaluating the law on real curve points):

    i1 = x_R   i2 = x_L   i3 = y_L   i4 = y_R      (L = left child, R = right child)
    i5 = x_out - Q        i6 = y_out

    N1 = 0  <=>  x_L + x_R + x_out = lambda^2  with lambda = B/A = (y_R-y_L)/(x_R-x_L)
    N2 = 0  <=>  y_L + y_out = lambda*(x_L - x_out)
i.e. exactly the chord addition law.  Q is a global additive offset on the x-slot;
consistency of the ladder up the tree forces Q = 0, so we take Q = 0 (and
ad_validate.py verifies that Q != 0 breaks the ladder, so this is not a free choice).

------------------------------------------------------------------------------
THE INTEGER-LIFT DISCIPLINE (this is the part the whole experiment is about)

u_k is a free integer, gcd(a_k,P) = gcd(c_k,P) = 1, so the exact integer content of
the congruence atom is

        c_k * P  |  a_k * L * Z_k          Z_k = c_k1*N1 + c_k2*N2
   <=>  P | L*Z_k            (the mod-P law: ordinary curve arithmetic)
   AND  c_k | a_k * L * Z_k  (the SMALL-MODULUS condition: invisible mod P)

On a live block the mod-P law holds, so Z_k = P * M_k with M_k an integer, and the
second condition is  m_k | M_k  where m_k = c_k / gcd(a_k, c_k).

M_k is NOT determined by the curve arithmetic: every slot is an integer variable
whose value is only fixed modulo P, so slot = alpha + P*t with alpha in [0,P) the
true reduced coordinate and t in Z a free "lift".  Writing

    i1 = a1+P*u1  i2 = a2+P*u2  i3 = a3+P*u3  i4 = a4+P*u4  i5 = a5+P*u5  i6 = a6+P*u6
    a = a1-a2  b = a4-a3  e = a1+a2+a5  g = a3+a6  h = a2-a5
    n1 = (e*a^2 - b^2)/P    n2 = (a*g - b*h)/P        (both integers -- checked)
    d = u1-u2   B_ = u4-u3   eps = u1+u2+u5   gam = u3+u6   eta = u2-u5

    nu1 := N1/P = n1 + (2*a*e*d + eps*a^2 - 2*b*B_) + P*(e*d^2 + 2*a*eps*d - B_^2)
                     + P^2*(eps*d^2)
    nu2 := N2/P = n2 + (a*gam + d*g - b*eta - B_*h) + P*(d*gam - B_*eta)

and the condition is  m_k | c_k1*nu1 + c_k2*nu2  for k = 1,2,3.

Differentiating gives  d(nu1,nu2)/d(u5,u6) = [[A^2,0],[B,A]] with A = a+P*d,
B = b+P*B_ -- exactly W's Jacobian, which is an independent check on the algebra:

    nu1 = nu1a + u5*A^2 ,   nu2 = nu2a + u5*B + u6*A .

Lift variables couple along the tree: a node's output lift IS its parent's input
lift.  So closure of the whole integer-lift family is a tree DP over lift residues.
By CRT the DP factors over the prime powers dividing the moduli, and each factor
is a DP over (Z/q)^2 states.  That is exact, not a heuristic.

Off-pins: on a dead block (L=0) the condition is c'_j*P | a'_j*i_j, i.e. P | i_j
(true value 0, the identity slot) and m'_j | (i_j/P).  The block's private output
slots i5,i6 are consumed ONLY through the mux multiplier cC = a*b = L, which is 0
exactly when the off-pin fires -- so i_j/P is a variable that occurs in no other
condition and w=0 always satisfies it.  See ad_validate.py, which checks this
structurally rather than asserting it.
"""
import os
from math import gcd

Q_OFFSET = 0  # forced; see module docstring


# ----------------------------------------------------------------- curve ----
class Curve:
    def __init__(self, p, a, b, N, G):
        self.p, self.a, self.b, self.N, self.G = p, a, b, N, tuple(G)

    def on(self, P):
        if P is None:
            return True
        x, y = P
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0

    def add(self, P1, P2):
        p, a = self.p, self.a
        if P1 is None:
            return P2
        if P2 is None:
            return P1
        x1, y1 = P1
        x2, y2 = P2
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1 + a) * pow(2 * y1 % p, -1, p) % p
        else:
            lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, k, P):
        R, Qp = None, P
        k %= self.N
        while k:
            if k & 1:
                R = self.add(R, Qp)
            Qp = self.add(Qp, Qp)
            k >>= 1
        return R


# ------------------------------------------------------------------ tree ----
class Node:
    __slots__ = ('idx', 'leaves', 'left', 'right', 'nid')

    def __init__(self, nid, idx=None, left=None, right=None):
        self.nid = nid
        self.idx = idx
        self.left = left
        self.right = right
        self.leaves = (idx,) if idx is not None else left.leaves + right.leaves


def build_tree(indices, mode='balanced', skew=0.6953125, counter=None):
    """Binary merge tree over the leaf indices.

    Shape:
      mode 'balanced'/'contig' -> split in half (real instance: depth 9 / 256 leaves)
      mode 'skew'              -> split 178:78 proportionally, the real root split.

    Leaf-index assignment: 'balanced' and 'skew' split the index set GREEDILY BY
    SUM (largest 2^i first into the lighter bucket), keeping the leaf counts at
    the shape's sizes.  This keeps every proper subtree's maximum partial sum
    below N, which is what makes the encoding faithful: a subtree sum reaching
    +-N mod N would make two merged children equal or opposite as curve points,
    which the gadget law cannot express (see CANCEL / DEGEN).  The real instance
    has the same property for free because N/2^256 = 1 - 2^-128.
    'contig' keeps the naive contiguous assignment, for comparison.
    """
    if counter is None:
        counter = [0]
    nid = counter[0]
    counter[0] += 1
    if len(indices) == 1:
        return Node(nid, idx=indices[0])
    if mode == 'skew':
        m = max(1, min(len(indices) - 1, int(round(len(indices) * skew))))
    else:
        m = len(indices) // 2
    if mode == 'contig':
        li, ri = indices[:m], indices[m:]
    else:
        li, ri, sl, sr = [], [], 0, 0
        for i in sorted(indices, reverse=True):
            if len(li) < m and (len(ri) >= len(indices) - m or sl <= sr):
                li.append(i)
                sl += 1 << i
            else:
                ri.append(i)
                sr += 1 << i
    l = build_tree(li, mode, skew, counter)
    r = build_tree(ri, mode, skew, counter)
    return Node(nid, left=l, right=r)


def postorder(root):
    out = []

    def rec(v):
        if v.idx is None:
            rec(v.left)
            rec(v.right)
        out.append(v)
    rec(root)
    return out


# ------------------------------------------------------- instance / eval ----
CHORD, DEGEN, CANCEL, PASS_L, PASS_R, DEAD = 0, 1, 2, 3, 4, 5
CLSNAME = {CHORD: 'chord', DEGEN: 'degeneracy(A=B=0)', CANCEL: 'cancel(A=0,B!=0)',
           PASS_L: 'pass-left', PASS_R: 'pass-right', DEAD: 'dead'}


class Instance:
    def __init__(self, curve, n, tree_mode='balanced'):
        self.cv = curve
        self.n = n
        self.ladder = []
        Ppt = curve.G
        for i in range(n):
            self.ladder.append(Ppt)
            Ppt = curve.add(Ppt, Ppt)
        self.root = build_tree(list(range(n)), tree_mode)
        self.nodes = postorder(self.root)
        self.internal = [v for v in self.nodes if v.idx is None]
        self.nnode = len(self.nodes)

    def k_of(self, S):
        return sum(1 << i for i in range(self.n) if (S >> i) & 1)

    def eval_subset(self, S):
        """Return (info per node, list of (node, cls, data)) for selector mask S.

        info[nid] = (live, k, point)  where point is the reduced (x,y) or None.
        """
        cv = self.cv
        info = {}
        blocks = []
        for v in self.nodes:
            if v.idx is not None:
                on = (S >> v.idx) & 1
                info[v.nid] = (bool(on), (1 << v.idx) if on else 0,
                               self.ladder[v.idx] if on else None)
                continue
            ll, kl, Pl = info[v.left.nid]
            lr, kr, Pr = info[v.right.nid]
            if ll and lr:
                if Pl[0] == Pr[0]:
                    cls = DEGEN if Pl[1] == Pr[1] else CANCEL
                else:
                    cls = CHORD
                Pv = cv.add(Pl, Pr)
                info[v.nid] = (True, kl + kr, Pv)
                blocks.append((v, cls, (Pl, Pr, Pv)))
            elif ll:
                info[v.nid] = (True, kl, Pl)
                blocks.append((v, PASS_L, None))
            elif lr:
                info[v.nid] = (True, kr, Pr)
                blocks.append((v, PASS_R, None))
            else:
                info[v.nid] = (False, 0, None)
                blocks.append((v, DEAD, None))
        return info, blocks


# --------------------------------------------------- coefficient drawing ----
class Coeffs:
    """Per-block coefficients of the law.

    C = [[c11,c12],[c21,c22],[c31,c32]] : the 3x2 integer matrix, rank 2, all
        three 2x2 minors nonzero (W verified this exhaustively on all 383 blocks).
    m = (m1,m2,m3) : the effective small moduli  m_k = c_k/gcd(a_k,c_k) of the
        integer-lift condition.  In the real instance 927 of 3707 lift conditions
        have c != 1 (25%), so m_k = 1 with prob 1-frac and otherwise a draw from
        the modulus pool.
    mo = (mo5,mo6) : same for the two off-pin conditions of the block.
    """

    def __init__(self, rng, nblocks, frac=0.25, pool=(2, 3, 4, 5, 7, 8, 9, 11, 13),
                 cmax=12):
        self.C, self.m, self.mo = [], [], []
        for _ in range(nblocks):
            while True:
                M = [[rng.randrange(-cmax, cmax + 1) for _ in range(2)] for _ in range(3)]
                mins = [M[0][0] * M[1][1] - M[0][1] * M[1][0],
                        M[0][0] * M[2][1] - M[0][1] * M[2][0],
                        M[1][0] * M[2][1] - M[1][1] * M[2][0]]
                if all(x != 0 for x in mins) and all(r[0] or r[1] for r in M):
                    break
            self.C.append(M)
            self.m.append(tuple(rng.choice(pool) if rng.random() < frac else 1
                                for _ in range(3)))
            self.mo.append(tuple(rng.choice(pool) if rng.random() < frac else 1
                                 for _ in range(2)))

    def prime_powers(self):
        """{prime: max exponent} over all congruence moduli."""
        out = {}
        for trip in self.m:
            for mm in trip:
                x, d = mm, 2
                while x > 1:
                    while d * d <= x and x % d:
                        d += 1
                    if d * d > x:
                        d = x
                    e = 0
                    while x % d == 0:
                        x //= d
                        e += 1
                    out[d] = max(out.get(d, 0), e)
                    d = 2
        return out


# ------------------------------------------------------------- the lift DP --
FULL_CACHE_LIMIT = 400000


class LiftDP:
    """Closure of the integer-lift (congruence) family, exactly, by tree DP.

    For one prime ell with max exponent E: q = ell**E, states (t,s) in (Z/q)^2
    encoded t*q+s, sets as Python int bitmasks.
    """

    def __init__(self, inst, co, ell, E, leaf_free=False, leaf_mod=None):
        self.inst, self.co = inst, co
        self.ell, self.E = ell, E
        self.q = ell ** E
        self.P = inst.cv.p % self.q
        self.full = (1 << (self.q * self.q)) - 1
        # per-block condition data reduced to this prime
        self.blk = {}
        for bi, v in enumerate(inst.internal):
            C = co.C[bi]
            mods = []
            for k in range(3):
                m = co.m[bi][k]
                e = 0
                while m % ell == 0:
                    m //= ell
                    e += 1
                mods.append(ell ** e)
            self.blk[v.nid] = (C, tuple(mods), any(x > 1 for x in mods))
        self.mask_cache = {}
        self.node_memo = [dict() for _ in range(inst.nnode)]
        self.leaf_free, self.leaf_mod = leaf_free, leaf_mod

    # --- leaf state
    def leaf_mask(self, idx, on):
        q = self.q
        if not self.leaf_free:
            return 1  # only (0,0)
        mt, ms = self.leaf_mod[idx]
        et = 1
        while mt % self.ell == 0:
            mt //= self.ell
            et *= self.ell
        es = 1
        ms2 = ms
        while ms2 % self.ell == 0:
            ms2 //= self.ell
            es *= self.ell
        msk = 0
        for t in range(0, q, et):
            for s in range(0, q, es):
                msk |= 1 << (t * q + s)
        return msk

    # --- the per-block solution mask, cached on (nid, A, B, nu1a, nu2a)
    def sol_mask(self, nid, A, B, nu1a, nu2a):
        key = (nid, A, B, nu1a, nu2a)
        m = self.mask_cache.get(key)
        if m is not None:
            return m
        q = self.q
        C, mods, _ = self.blk[nid]
        A2 = A * A % q
        msk = 0
        for u5 in range(q):
            n1 = (nu1a + u5 * A2) % q
            n2b = (nu2a + u5 * B) % q
            for u6 in range(q):
                n2 = (n2b + u6 * A) % q
                ok = True
                for k in range(3):
                    mk = mods[k]
                    if mk == 1:
                        continue
                    if (C[k][0] * n1 + C[k][1] * n2) % mk:
                        ok = False
                        break
                if ok:
                    msk |= 1 << (u5 * q + u6)
        if len(self.mask_cache) < FULL_CACHE_LIMIT:
            self.mask_cache[key] = msk
        return msk

    def states(self, mask):
        q = self.q
        out = []
        i = 0
        while mask:
            if mask & 1:
                out.append((i // q, i % q))
            mask >>= 1
            i += 1
        return out

    def run(self, S, info=None):
        """Return the root feasible-lift bitmask (0 == the family does NOT close)."""
        inst = self.inst
        if info is None:
            info, _ = inst.eval_subset(S)
        q, P = self.q, self.P
        cur = {}
        for v in inst.nodes:
            if v.idx is not None:
                cur[v.nid] = self.leaf_mask(v.idx, (S >> v.idx) & 1)
                continue
            pat = 0
            for j, ix in enumerate(v.leaves):
                if (S >> ix) & 1:
                    pat |= 1 << j
            memo = self.node_memo[v.nid]
            got = memo.get(pat)
            if got is not None:
                cur[v.nid] = got
                continue
            ll, _, Pl = info[v.left.nid]
            lr, _, Pr = info[v.right.nid]
            if not ll and not lr:
                res = 1
            elif ll and not lr:
                res = cur[v.left.nid]
            elif lr and not ll:
                res = cur[v.right.nid]
            else:
                if Pl[0] == Pr[0]:
                    # DEGEN: law vacuous, output completely free.
                    # CANCEL: block is integrally infeasible -> empty.
                    res = self.full if Pl[1] == Pr[1] else 0
                else:
                    Pv = inst.cv.add(Pl, Pr)
                    res = self.trans(v.nid, Pl, Pr, Pv,
                                     cur[v.left.nid], cur[v.right.nid])
            memo[pat] = res
            cur[v.nid] = res
            if res == 0:
                # empty at a node => empty at the root (tree DP, no way back)
                return 0
        return cur[inst.root.nid]

    def trans(self, nid, Pl, Pr, Pv, Fl, Fr):
        q, P = self.q, self.P
        Pn = self.inst.cv.p
        xl, yl = Pl
        xr, yr = Pr
        xv, yv = Pv
        aa = xr - xl
        bb = yr - yl
        ee = xr + xl + xv + Q_OFFSET
        gg = yl + yv
        hh = xl - xv
        num1 = ee * aa * aa - bb * bb
        num2 = aa * gg - bb * hh
        assert num1 % Pn == 0 and num2 % Pn == 0, 'law not satisfied mod P'
        n1 = (num1 // Pn) % q
        n2 = (num2 // Pn) % q
        aa %= q
        bb %= q
        ee %= q
        gg %= q
        hh %= q
        Ls = self.states(Fl)
        Rs = self.states(Fr)
        res = 0
        full = self.full
        sol = self.sol_mask
        for (tl, sl) in Ls:
            for (tr, sr) in Rs:
                d = (tr - tl) % q
                bt = (sr - sl) % q
                eps0 = (tl + tr) % q
                A = (aa + P * d) % q
                B = (bb + P * bt) % q
                nu1a = (n1 + 2 * aa * ee * d + eps0 * aa * aa - 2 * bb * bt
                        + P * (ee * d * d + 2 * aa * eps0 * d - bt * bt)
                        + P * P * eps0 * d * d) % q
                nu2a = (n2 + aa * sl + d * gg - bb * tl - bt * hh
                        + P * (d * sl - bt * tl)) % q
                res |= sol(nid, A, B, nu1a, nu2a)
                if res == full:
                    return res
        return res


# --------------------------------------------------------------- closure ----
def congruence_closure(inst, co, leaf_free=False, leaf_mod=None, dps=None):
    """Build (or reuse) one LiftDP per prime dividing the moduli."""
    if dps is None:
        pp = co.prime_powers()
        dps = [LiftDP(inst, co, ell, E, leaf_free, leaf_mod) for ell, E in sorted(pp.items())]
    return dps


def closes(dps, S, info=None):
    for dp in dps:
        if dp.run(S, info) == 0:
            return False
    return True


def root_pinned_closes(dps, S, info=None):
    for dp in dps:
        m = dp.run(S, info)
        if not (m & 1):
            return False
    return True
