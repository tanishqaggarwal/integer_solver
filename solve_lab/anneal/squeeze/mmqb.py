#!/usr/bin/env python3
"""mmqb.py -- QUBO compiler extensions for the modmul squeeze.

Everything here is additive on top of ../qubo.py:QB.  Three things the base
compiler cannot express, each of which turns out to matter:

  1. `assert_terms` -- assert an integer identity given as a list of SIGNED
     POWER-OF-TWO terms  (monomial, sign, shift)  instead of a dict
     monomial -> integer coefficient.  The base `assert_zero` binary-expands
     each coefficient, so a coefficient of  p = 2^256-2^32-977  contributes
     popcount(p) = 250 column entries.  Written as a signed sum it is 6.
     This is the whole pseudo-Mersenne saving, and it generalises: `naf_split`
     finds the shortest signed representation of any constant.

  2. Extra carry disciplines: 'dadda'  (3:2 tree stopped at height 2 with a
     final ripple), 'unary' (thermometer-coded carries).

  3. A word type that carries its own witness evaluator, so Karatsuba /
     Toom-Cook trees can be replayed by QB.witness unchanged.
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qubo import QB                                              # noqa: E402


# ----------------------------------------------------------- constant splits
def bin_split(n):
    """plain binary expansion of |n| with the sign of n: [(sign, shift), ...]"""
    s = 1 if n >= 0 else -1
    n = abs(n)
    return [(s, t) for t in range(n.bit_length()) if (n >> t) & 1]


def naf_split(n):
    """non-adjacent form: the shortest signed-binary representation.
       For p = 2^256-2^32-977 this returns 6 terms instead of 250."""
    s = 1 if n >= 0 else -1
    n = abs(n)
    out, t = [], 0
    while n:
        if n & 1:
            z = 2 - (n % 4)
            out.append((s if z > 0 else -s, t))
            n -= z
        n >>= 1
        t += 1
    return out


def best_split(n):
    a, b = bin_split(n), naf_split(n)
    return a if len(a) <= len(b) else b


# ------------------------------------------------------------------- words
class Wrd:
    """a bit vector plus a witness evaluator wv -> int."""
    __slots__ = ('bits', 'val', 'name')

    def __init__(self, bits, val, name=''):
        self.bits, self.val, self.name = bits, val, name

    def __len__(self):
        return len(self.bits)

    def slice(self, lo, hi):
        b = self.bits[lo:hi]
        f = self.val
        return Wrd(b, lambda wv, f=f, lo=lo, m=(1 << (hi - lo)) - 1: (f(wv) >> lo) & m)


class MMQB(QB):
    """QB + signed-term identities + extra carry disciplines."""

    def __init__(self, chunk=16, mode='wallace', dadda_height=2):
        super().__init__(chunk=chunk, mode=mode)
        self.dadda_height = dadda_height
        self.and_lookups = 0        # AND() calls
        self.and_hits = 0           # AND() calls served from the cache
        self.squares = []           # every (lin, const) asserted as a square == 0
        self.orders = []            # (u, v) meaning v <= u  (thermometer carries)

    # --- record the equations, for the presolver -------------------------
    def add_square(self, lin, const, tgt=None):
        if tgt is None:
            self.squares.append((dict(lin), const))
        return super().add_square(lin, const, tgt)

    # --- instrumented AND cache -----------------------------------------
    def AND(self, i, j):
        self.and_lookups += 1
        if i == j:
            self.and_hits += 1
            return i
        key = (i, j) if i < j else (j, i)
        if key in self.andcache:
            self.and_hits += 1
            return self.andcache[key]
        return super().AND(i, j)

    # --- a word that knows its own value ---------------------------------
    def mkword(self, name, nbits, fn):
        bits = self.word(name, nbits, fn)
        return Wrd(bits, (lambda wv, name=name: wv[name]), name)

    # --- the signed-term identity ----------------------------------------
    def assert_terms(self, terms, consts, tag):
        """assert  sum_k sign_k * 2^shift_k * mono_k  +  sum consts  == 0  over Z.
           terms:  iterable of (monomial tuple, sign, shift)
           consts: iterable of (sign, shift)"""
        cols_pos, cols_neg = defaultdict(list), defaultdict(list)
        for mono, sg, sh in terms:
            v = self.mono_var(mono)
            (cols_pos if sg > 0 else cols_neg)[sh].append(v)
        kpos = kneg = 0
        for sg, sh in consts:
            if sg > 0:
                kpos += 1 << sh
            else:
                kneg += 1 << sh
        d = min(kpos, kneg)
        kpos -= d
        kneg -= d
        if self.mode == 'wallace':
            return self._wallace_eq(cols_pos, kpos, cols_neg, kneg, tag)
        if self.mode == 'dadda':
            return self._dadda_eq(cols_pos, kpos, cols_neg, kneg, tag)
        if self.mode == 'unary':
            return self._unary_eq(cols_pos, kpos, cols_neg, kneg, tag)
        return self._ripple_eq(cols_pos, kpos, cols_neg, kneg, tag)

    # --- mode C: 3:2 tree down to height H, then a binary ripple at the top --
    def _reduce_to(self, cols, kconst, tag, height):
        work = defaultdict(list)
        for c, l in cols.items():
            work[c] += list(l)
        t = 0
        while kconst >> t:
            if (kconst >> t) & 1:
                work[t].append(self.ONE)
            t += 1
        c, top = 0, (max(work) if work else 0)
        while c <= top:
            while len(work[c]) > height:
                grp = [work[c].pop() for _ in range(min(3, len(work[c])))]
                nc = sum(1 for g in grp if g == self.ONE)
                vs = [g for g in grp if g != self.ONE]
                if len(grp) == 2 and nc == 2:
                    work[c + 1].append(self.ONE)
                    top = max(top, c + 1)
                    continue
                sv = self.new(f"fa:{tag}:{c}:s", 'adder')
                dv = self.new(f"fa:{tag}:{c}:c", 'adder')
                self.n_carry += 2
                lin = defaultdict(int)
                for v in vs:
                    lin[v] += 1
                lin[sv] -= 1
                lin[dv] -= 2
                self.add_square(dict(lin), nc)
                self.trace.append(('fa', vs, nc, sv, dv))
                work[c].append(sv)
                work[c + 1].append(dv)
                top = max(top, c + 1)
            c += 1
            top = max(top, max(work) if work else 0)
        return work

    def _dadda_eq(self, cols_pos, kpos, cols_neg, kneg, tag):
        H = self.dadda_height
        A = self._reduce_to(cols_pos, kpos, tag + ':P', H)
        B = self._reduce_to(cols_neg, kneg, tag + ':N', H)
        # final ripple: one balanced binary-carry equation per column
        cp, cn = defaultdict(list), defaultdict(list)
        kp = kn = 0
        for c in set(A) | set(B):
            for v in A.get(c, []):
                if v == self.ONE:
                    kp += 1 << c
                else:
                    cp[c].append(v)
            for v in B.get(c, []):
                if v == self.ONE:
                    kn += 1 << c
                else:
                    cn[c].append(v)
        self._ripple_eq(cp, kp, cn, kn, tag + ':R')

    # --- mode D: unary / thermometer carries ------------------------------
    def _unary_eq(self, cols_pos, kpos, cols_neg, kneg, tag):
        ncol = max(list(cols_pos) + list(cols_neg) +
                   [kpos.bit_length(), kneg.bit_length()] + [0]) + 1
        rec = []
        cin_vars, cin_base = [], 0
        for c in range(ncol + 4):
            pos, neg = cols_pos.get(c, []), cols_neg.get(c, [])
            kp, kn = (kpos >> c) & 1, (kneg >> c) & 1
            if not pos and not neg and not kp and not kn and not cin_vars and cin_base == 0:
                break
            lin = defaultdict(int)
            for v in pos:
                lin[v] += 1
            for v in neg:
                lin[v] -= 1
            for v in cin_vars:
                lin[v] += 1
            const_c = kp - kn + cin_base
            lo = sum(min(0, w) for w in lin.values()) + const_c
            hi = sum(max(0, w) for w in lin.values()) + const_c
            co_lo, co_hi = -((-lo) // 2), hi // 2
            if co_lo > co_hi:
                co_lo = co_hi = 0
            nb = co_hi - co_lo                 # unary: one variable per level
            cvars = [self.new(f"ucarry:{tag}:{c}:{t}", 'carry') for t in range(nb)]
            self.n_carry += nb
            for t in range(nb - 1):            # thermometer order: c_{t+1} <= c_t
                self.pen[(cvars[t + 1],)] += 1
                self.pen[(min(cvars[t], cvars[t + 1]), max(cvars[t], cvars[t + 1]))] += -1
                self.orders.append((cvars[t], cvars[t + 1]))
                self.max_clique = max(self.max_clique, 2)
            lin2 = dict(lin)
            for v in cvars:
                lin2[v] = lin2.get(v, 0) - 2
            self.add_square(lin2, const_c - 2 * co_lo)
            rec.append((dict(lin), const_c, cvars, co_lo, 'unary'))
            cin_vars = list(cvars)
            cin_base = co_lo
        self.trace.append(('ucols', rec))

    # --- witness support for unary columns -------------------------------
    def witness(self, inputs, wv0=None):
        x = [0] * self.n
        for v, val in inputs.items():
            x[v] = val
        wv = dict(wv0 or {})
        for op in self.trace:
            if op[0] == 'ucols':
                for lin, const_c, cvars, base, _ in op[1]:
                    S = const_c + sum(w * x[v] for v, w in lin.items())
                    assert S % 2 == 0, "column parity violated"
                    co = S // 2
                    assert base <= co <= base + len(cvars), "carry out of range"
                    for t in range(len(cvars)):
                        x[cvars[t]] = 1 if t < co - base else 0
                continue
            if op[0] == 'and':
                _, z, i, j = op
                x[z] = x[i] * x[j]
            elif op[0] == 'word':
                _, name, bits, fn = op
                val = fn(wv)
                wv[name] = val
                for t, v in enumerate(bits):
                    x[v] = (val >> t) & 1
                assert val >= 0 and val >> len(bits) == 0, f"word {name} overflow ({val})"
            elif op[0] == 'fa':
                _, vs, nc, sv, dv = op
                tot = nc + sum(x[v] for v in vs)
                x[sv] = tot & 1
                x[dv] = tot >> 1
            elif op[0] == 'chunk':
                _, chunk, avars = op
                s = sum(x[v] for v in chunk)
                for t, v in enumerate(avars):
                    x[v] = (s >> t) & 1
            elif op[0] == 'cols':
                for lin, const_c, cvars, base in op[1]:
                    S = const_c + sum(w * x[v] for v, w in lin.items())
                    assert S % 2 == 0, "column parity violated"
                    co = S // 2
                    assert co - base >= 0 and (co - base) >> len(cvars) == 0, "carry out of range"
                    for t, v in enumerate(cvars):
                        x[v] = ((co - base) >> t) & 1
        return x, wv

    # --- soundness of the AND weight -------------------------------------
    def and_weight_ok(self):
        """finalize() picks W = 1 + max over AND vars of their |coefficient| load
        in the non-AND part.  Re-check it here: flipping an AND var away from
        its product changes the AND penalty by >= W and the rest by < W, so
        every E=0 state has every AND var correct."""
        andvars = set(self.andcache.values())
        load = defaultdict(int)
        for m, c in self.pen.items():
            for v in m:
                if v in andvars:
                    load[v] += abs(c)
        return (max(load.values()) if load else 0) < self.W

    def triple(self):
        st = self.stats()
        return st['vars'], st['max_clique'], st['dynamic_range_bits']
