"""EXACT integer polynomial system for the witness region — products carried, nothing probed away.

Everything measured on this thread so far went through finite differences, which can only ever see
an affine model.  Here the knobs are replaced by honest symbols `t_j` and the DAG is forward-
evaluated with exact multivariate polynomial arithmetic over Z, so every `x_a * x_b` in an atom is
carried as a product.  The 12 region rows come out as exact elements of Z[t_1..t_k]; no square is
stripped, no difference is taken.

Usage:  python3 polyexact.py [narrow|wide]
"""
import os, sys, json, time, itertools
from collections import defaultdict

import ev
import optN
from optN import make, build, WIT, POOL, fr, FREE, FR0, atom_eqs, _bits
import frameB as FB

HERE = os.path.dirname(os.path.abspath(__file__))
eq_terms = ev.eq_terms


# ---------------------------------------------------------------- exact sparse polynomial over Z
class P:
    """Multivariate polynomial over Z in NK variables; monomials keyed by exponent tuple."""
    __slots__ = ('c',)
    NK = 0

    def __init__(self, c=None):
        self.c = c or {}

    @staticmethod
    def const(n):
        return P({(0,) * P.NK: n} if n else {})

    @staticmethod
    def var(j, base=0):
        e = [0] * P.NK
        e[j] = 1
        d = {tuple(e): 1}
        if base:
            d[(0,) * P.NK] = base
        return P(d)

    def __bool__(self):
        return bool(self.c)

    def __eq__(self, o):
        if isinstance(o, int):
            return self.c == ({(0,) * P.NK: o} if o else {})
        return self.c == o.c

    def _add(self, o, sgn):
        r = dict(self.c)
        for k, v in o.c.items():
            n = r.get(k, 0) + sgn * v
            if n:
                r[k] = n
            elif k in r:
                del r[k]
        return P(r)

    def __add__(self, o):
        return self._add(o if isinstance(o, P) else P.const(o), 1)

    def __radd__(self, o):
        return self.__add__(o)

    def __sub__(self, o):
        return self._add(o if isinstance(o, P) else P.const(o), -1)

    def __rsub__(self, o):
        return (P.const(o) if isinstance(o, int) else o)._add(self, -1)

    def __neg__(self):
        return P({k: -v for k, v in self.c.items()})

    def __mul__(self, o):
        if isinstance(o, int):
            return P({k: v * o for k, v in self.c.items()}) if o else P()
        r = {}
        for k1, v1 in self.c.items():
            for k2, v2 in o.c.items():
                k = tuple(a + b for a, b in zip(k1, k2))
                n = r.get(k, 0) + v1 * v2
                if n:
                    r[k] = n
                elif k in r:
                    del r[k]
        return P(r)

    def __rmul__(self, o):
        return self.__mul__(o)

    def __pow__(self, n):
        r = P.const(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def deg(self):
        return max((sum(k) for k in self.c), default=-1)

    def nterms(self):
        return len(self.c)

    def maxbits(self):
        return max((abs(v).bit_length() for v in self.c.values()), default=0)

    def support_vars(self):
        s = set()
        for k in self.c:
            for j, e in enumerate(k):
                if e:
                    s.add(j)
        return s

    def const_term(self):
        return self.c.get((0,) * P.NK, 0)


# ---------------------------------------------------------------- build the state and the knobs
def region_and_knobs(mode):
    st = make(WIT)
    b0 = build(st)
    Rl = b0['R']
    if mode == 'narrow':
        knobs = b0['knobs']
    else:
        import widen
        knobs, outside = widen.wide_knobs(st, Rl, verbose=True)
    return st, Rl, knobs


def symbolic_rows(st, Rl, knobs, verbose=True):
    """Exact polynomial value of each region row's inner sum, in Z[t_0..t_{k-1}]."""
    P.NK = len(knobs)
    v = list(st.v)                      # numeric everywhere
    ns = {'v': v, '__builtins__': {}}

    # place symbols on the knobs
    aff = set()
    ck = set()
    for j, Y in enumerate(knobs):
        v[Y] = P.var(j, st.fv.get(Y, 0))
        aff.update(fr.desc[Y])
        ck.update(fr.chk[Y])
    if verbose:
        print('knobs: %d   downstream DAG variables touched: %d   check atoms touched: %d'
              % (len(knobs), len(aff), len(ck)), flush=True)

    # forward-evaluate the affected definitions in topological order, exactly
    t0 = time.time()
    symvars = []
    for u in sorted(aff, key=lambda u: fr.pos[u]):
        v[u] = eval(FB.DEFEXPR[u], ns)
        if isinstance(v[u], P):
            symvars.append(u)
    if verbose:
        print('   forward pass %.1fs ; variables that are genuinely symbolic: %d'
              % (time.time() - t0, len(symvars) + len(knobs)), flush=True)
        if symvars:
            dd = [v[u].deg() for u in symvars]
            print('   downstream variable degrees: min %d max %d' % (min(dd), max(dd)), flush=True)

    # atom values
    av = dict(st.av)
    for a in sorted(ck):
        av[a] = eval(FB.ACODE[a], ns)

    # rows: the equation inner sum, EXACT (no square stripping, no linear core)
    rows = {}
    for e in Rl:
        m, sq, tl = eq_terms[e]
        acc = P()
        for c, a in tl:
            x = av.get(a)
            if isinstance(x, P):
                acc = acc + x * c
            elif x:
                acc = acc + P.const(c * x)
        rows[e] = (acc, sq, m)
    return rows, av, symvars, v


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'narrow'
    st, Rl, knobs = region_and_knobs(mode)
    print('=== EXACT POLYNOMIAL SIZING (%s knob set) ===' % mode, flush=True)
    print('region |R| = %d : %s' % (len(Rl), Rl), flush=True)
    print('knobs (%d): %s' % (len(knobs), knobs), flush=True)

    rows, av, symvars, v = symbolic_rows(st, Rl, knobs)

    print('\nper-row exact polynomial (inner sum; `sq` = the equation squares it):', flush=True)
    out = []
    maxdeg = 0
    for e in Rl:
        pol, sq, m = rows[e]
        d = pol.deg()
        maxdeg = max(maxdeg, d if not sq else 2 * d)
        supp = sorted(pol.support_vars())
        print('   eq %-6d  deg %-3d  terms %-6d  maxcoefbits %-6d  sq=%-5s  const=%s  knobs used %s'
              % (e, d, pol.nterms(), pol.maxbits(), sq,
                 ('0' if pol.const_term() == 0 else '%d-bit' % pol.const_term().bit_length()),
                 [knobs[j] for j in supp]), flush=True)
        out.append(dict(eq=e, deg=d, terms=pol.nterms(), bits=pol.maxbits(), sq=bool(sq),
                        const_bits=(0 if pol.const_term() == 0
                                    else abs(pol.const_term()).bit_length()),
                        knobs_used=[knobs[j] for j in supp]))

    print('\n=== system size ===', flush=True)
    print('unknowns %d, equations %d, max total degree %d' % (len(knobs), len(Rl), maxdeg),
          flush=True)
    json.dump(dict(mode=mode, knobs=knobs, R=Rl, rows=out, maxdeg=maxdeg,
                   symvars=len(symvars) + len(knobs)),
              open(os.path.join(HERE, 'runs', 'polyexact_%s.json' % mode), 'w'), indent=1)
    return rows, knobs, Rl, st


if __name__ == '__main__':
    main()
