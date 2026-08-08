#!/usr/bin/env python3
"""
qubo.py -- a QUBO/Ising compiler for the decision core of EQUATIONS.txt.

The core (see structure.py) is:  find b_0..b_{m-1} in {0,1} with
    sum_i b_i * P_i  ==  T      on  E : y^2 = x^3 + B  over F_p,   P_i = 2^i G.

Encoding strategy (see ENCODING.md):
  * signed-digit rewrite  e_i = 2 b_i - 1  so every ladder step is an
    UNCONDITIONAL addition of the *constant* point P_i with a sign that is
    linear in b_i.  x-coordinate of the addend is then a compile-time constant
    and only its y-coordinate carries a variable (a single b_i).
  * affine addition with an explicit inverse witness  d*dinv = 1, which both
    removes the fraction-free cubic AND closes the degenerate-division
    loophole (d = 0) for free -- 4 modular multiplications per step instead of
    5 fraction-free + 1 invertibility check.
  * every value is a normalised s-bit word; every relation is asserted as an
    exact integer identity  LHS - RHS - p*q = 0  with an explicit quotient
    word q, then balanced column-by-column with bounded carries.  No 2^{2s}
    coefficients ever appear in the Hamiltonian, so the coupler dynamic range
    stays polynomial instead of exponential.
  * tall columns are chunk-compressed (parameter K) to trade a few ancillas
    for a quadratic reduction in coupler count.

Every penalty is a sum of squares of integer-valued linear forms plus AND
penalties, so  E >= 0  always and  E == 0  exactly on solutions.
"""
from collections import defaultdict


# ----------------------------------------------------------------- QUBO core
class QB:
    def __init__(self, chunk=16, mode='binary'):
        self.mode = mode          # 'binary' (fewest qubits) | 'wallace' (lowest |J| range)
        self.n = 0
        self.names = []
        self.kind = []
        self.pen = defaultdict(int)      # monomial (deg<=2) -> coefficient
        self.andpen = defaultdict(int)   # AND penalties, unit weight
        self.andcache = {}
        self.trace = []                  # ops replayed by witness()
        self.chunk = chunk
        self.W_and = None
        self.max_clique = 0
        self.n_and = 0
        self.n_carry = 0
        self.n_word = 0

    def new(self, name, kind):
        self.names.append(name); self.kind.append(kind)
        self.n += 1
        return self.n - 1

    # --- penalty helpers -------------------------------------------------
    def _add(self, tgt, mono, c):
        if c: tgt[mono] += c

    def add_square(self, lin, const, tgt=None):
        """add (sum lin[v]*v + const)^2 ; v binary so v^2 = v.
        Each such square makes a clique on its support, so the largest one is the
        binding constraint for minor-embedding: a K_c needs ~c^2/4 physical qubits
        on a Pegasus/Zephyr graph."""
        tgt = self.pen if tgt is None else tgt
        if tgt is self.pen and len(lin) > self.max_clique: self.max_clique = len(lin)
        vs = sorted(lin)
        self._add(tgt, (), const * const)
        for v in vs:
            c = lin[v]
            self._add(tgt, (v,), c * c + 2 * const * c)
        for a in range(len(vs)):
            for b in range(a + 1, len(vs)):
                self._add(tgt, (vs[a], vs[b]), 2 * lin[vs[a]] * lin[vs[b]])

    def AND(self, i, j):
        if i == j: return i
        key = (i, j) if i < j else (j, i)
        if key in self.andcache: return self.andcache[key]
        z = self.new(f"and({i},{j})", 'and'); self.n_and += 1
        # a*b - 2az - 2bz + 3z  >= 0, == 0 iff z = a*b
        self.andpen[key] += 1
        self.andpen[(min(i, z), max(i, z))] += -2
        self.andpen[(min(j, z), max(j, z))] += -2
        self.andpen[(z,)] += 3
        self.andcache[key] = z
        self.trace.append(('and', z, i, j))
        return z

    def mono_var(self, mono):
        """collapse a monomial (tuple of var ids) to a single variable."""
        mono = tuple(sorted(set(mono)))
        if len(mono) == 0: return None
        v = mono[0]
        for u in mono[1:]:
            v = self.AND(v, u)
        return v

    # --- words -----------------------------------------------------------
    def word(self, name, nbits, fn):
        """fresh nbits-wide unsigned word; fn(wv) -> its integer value (witness)."""
        bits = [self.new(f"{name}[{t}]", 'word') for t in range(nbits)]
        self.n_word += nbits
        self.trace.append(('word', name, bits, fn))
        return bits

    # --- the workhorse: assert an integer polynomial identity == 0 --------
    def assert_zero(self, poly, const, tag):
        """poly: dict monomial(tuple of var ids) -> int coefficient.
           Enforces sum poly + const == 0 over Z."""
        cols_pos, cols_neg = defaultdict(list), defaultdict(list)
        for mono, c in poly.items():
            if c == 0: continue
            v = self.mono_var(mono)
            tgt = cols_pos if c > 0 else cols_neg
            a, t = abs(c), 0
            while a:
                if a & 1: tgt[t].append(v)
                a >>= 1; t += 1
        kpos = const if const > 0 else 0
        kneg = -const if const < 0 else 0
        if self.mode == 'wallace':
            return self._wallace_eq(cols_pos, kpos, cols_neg, kneg, tag)
        return self._ripple_eq(cols_pos, kpos, cols_neg, kneg, tag)

    # --- mode A: one balanced equation per column, carries in binary (min qubits) ---
    def _ripple_eq(self, cols_pos, kpos, cols_neg, kneg, tag):
        ncol = max(list(cols_pos) + list(cols_neg) + [kpos.bit_length(), kneg.bit_length()] + [0]) + 1
        rec = []
        cin_vars, cin_base, cin_lo, cin_hi = [], 0, 0, 0
        for c in range(ncol + 4):
            pos, neg = cols_pos.get(c, []), cols_neg.get(c, [])
            kp, kn = (kpos >> c) & 1, (kneg >> c) & 1
            if not pos and not neg and not kp and not kn and cin_lo == 0 and cin_hi == 0:
                break
            lin = defaultdict(int)
            for v, w in self._compress(pos, +1, tag, c): lin[v] += w
            for v, w in self._compress(neg, -1, tag, c): lin[v] += w
            for v, w in cin_vars: lin[v] += w
            const_c = kp - kn + cin_base
            lo = sum(min(0, w) for w in lin.values()) + const_c
            hi = sum(max(0, w) for w in lin.values()) + const_c
            co_lo, co_hi = -((-lo) // 2), hi // 2
            if co_lo > co_hi: co_lo = co_hi = 0
            nb = max(0, (co_hi - co_lo).bit_length())
            cvars = [self.new(f"carry:{tag}:{c}:{t}", 'carry') for t in range(nb)]
            self.n_carry += nb
            lin2 = dict(lin)
            for t, v in enumerate(cvars): lin2[v] = lin2.get(v, 0) - (1 << (t + 1))
            self.add_square(lin2, const_c - 2 * co_lo)
            rec.append((dict(lin), const_c, cvars, co_lo))
            cin_vars = [(v, (1 << t)) for t, v in enumerate(cvars)]
            cin_base, cin_lo = co_lo, co_lo
            cin_hi = co_lo + ((1 << nb) - 1 if nb else 0)
        self.trace.append(('cols', rec))

    # --- mode B: Wallace-tree compression, every penalty has |coef| <= 4 ---
    ONE = -1

    def _wallace(self, cols, kconst, tag):
        """compress a weighted bag to at most one term per column using 3:2 / 2:2 adders."""
        work = defaultdict(list)
        for c, l in cols.items(): work[c] += list(l)
        t = 0
        while kconst >> t:
            if (kconst >> t) & 1: work[t].append(self.ONE)
            t += 1
        c = 0
        top = max(work) if work else 0
        while c <= top:
            while len(work[c]) > 1:
                grp = [work[c].pop() for _ in range(min(3, len(work[c])))]
                nc = sum(1 for g in grp if g == self.ONE)
                vs = [g for g in grp if g != self.ONE]
                if len(grp) == 2 and nc == 2:
                    work[c + 1].append(self.ONE); top = max(top, c + 1); continue
                sv = self.new(f"fa:{tag}:{c}:s", 'adder')
                dv = self.new(f"fa:{tag}:{c}:c", 'adder')
                self.n_carry += 2
                lin = defaultdict(int)
                for v in vs: lin[v] += 1
                lin[sv] -= 1; lin[dv] -= 2
                self.add_square(dict(lin), nc)
                self.trace.append(('fa', vs, nc, sv, dv))
                work[c].append(sv); work[c + 1].append(dv)
                top = max(top, c + 1)
            c += 1
            top = max(top, max(work) if work else 0)
        return {c: (l[0] if l else None) for c, l in work.items()}

    def _wallace_eq(self, cols_pos, kpos, cols_neg, kneg, tag):
        A = self._wallace(cols_pos, kpos, tag + ':P')
        B = self._wallace(cols_neg, kneg, tag + ':N')
        for c in set(A) | set(B):
            a, b = A.get(c), B.get(c)
            lin, k = {}, 0
            for term, sg in ((a, 1), (b, -1)):
                if term is None: continue
                if term == self.ONE: k += sg
                else: lin[term] = lin.get(term, 0) + sg
            if not lin and k == 0: continue
            self.add_square(lin, k)

    def _compress(self, terms, sign, tag, col):
        """replace a long list of unit terms by binary-encoded chunk sums."""
        K = self.chunk
        if len(terms) <= K:
            return [(v, sign) for v in terms]
        out = []
        for g in range(0, len(terms), K):
            chunk = terms[g:g + K]
            if len(chunk) <= 2:
                out += [(v, sign) for v in chunk]; continue
            nb = len(chunk).bit_length()
            avars = [self.new(f"chunk:{tag}:{col}:{g}:{t}", 'chunk') for t in range(nb)]
            self.n_carry += nb
            lin = {v: 1 for v in chunk}
            for t, v in enumerate(avars): lin[v] = lin.get(v, 0) - (1 << t)
            self.add_square(lin, 0)
            self.trace.append(('chunk', chunk, avars))
            out += [(v, sign * (1 << t)) for t, v in enumerate(avars)]
        return out

    # --- assembly --------------------------------------------------------
    def finalize(self):
        """merge AND penalties with the smallest weight that provably cannot be traded away.

        Flipping one variable changes the non-AND part of the energy by at most
        the sum of |coefficients| of the penalty terms it occurs in; an AND
        penalty pays >= W for being wrong, so W above that local bound makes
        every AND gate rigid.  This is ~2^30 tighter than a global bound and is
        what keeps the coupler dynamic range inside an annealer's precision."""
        andvars = {z for _, z in self.andcache.items()}
        load = defaultdict(int)
        for m, c in self.pen.items():
            for v in m:
                if v in andvars: load[v] += abs(c)
        W = self.W_and or ((max(load.values()) if load else 0) + 1)
        Q = defaultdict(int)
        for m, c in self.pen.items(): Q[m] += c
        for m, c in self.andpen.items(): Q[m] += W * c
        self.W = W
        self.Q = {m: c for m, c in Q.items() if c}
        return self.Q

    def energy(self, x):
        e = 0
        for m, c in self.Q.items():
            if not m: e += c
            elif len(m) == 1: e += c * x[m[0]]
            else: e += c * x[m[0]] * x[m[1]]
        return e

    def stats(self):
        lin = sum(1 for m in self.Q if len(m) == 1)
        quad = sum(1 for m in self.Q if len(m) == 2)
        mags = [abs(c) for m, c in self.Q.items() if m]
        return dict(vars=self.n, linear=lin, couplers=quad, max_clique=self.max_clique,
                    and_vars=self.n_and, word_bits=self.n_word, carry_bits=self.n_carry,
                    max_coef=max(mags), min_coef=min(mags),
                    dynamic_range_bits=(max(mags) // max(1, min(mags))).bit_length())

    # --- witness ---------------------------------------------------------
    def witness(self, inputs, wv0=None):
        """replay the construction to fill every ancilla from the input bits."""
        x = [0] * self.n
        for v, val in inputs.items(): x[v] = val
        wv = dict(wv0 or {})
        for op in self.trace:
            if op[0] == 'and':
                _, z, i, j = op; x[z] = x[i] * x[j]
            elif op[0] == 'word':
                _, name, bits, fn = op
                val = fn(wv); wv[name] = val
                for t, v in enumerate(bits): x[v] = (val >> t) & 1
                assert val >> len(bits) == 0, f"word {name} overflow"
            elif op[0] == 'fa':
                _, vs, nc, sv, dv = op
                tot = nc + sum(x[v] for v in vs)
                x[sv] = tot & 1; x[dv] = tot >> 1
            elif op[0] == 'chunk':
                _, chunk, avars = op
                s = sum(x[v] for v in chunk)
                for t, v in enumerate(avars): x[v] = (s >> t) & 1
            elif op[0] == 'cols':
                for lin, const_c, cvars, base in op[1]:
                    S = const_c + sum(w * x[v] for v, w in lin.items())
                    assert S % 2 == 0, "column parity violated"
                    co = S // 2
                    assert co - base >= 0 and (co - base) >> len(cvars) == 0, "carry out of range"
                    for t, v in enumerate(cvars): x[v] = ((co - base) >> t) & 1
        return x, wv
