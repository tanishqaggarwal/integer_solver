#!/usr/bin/env python3
"""lp_core.py -- extract the linear+AND constraint system of a modmul QUBO and
build its LP relaxation, for convex-structure / persistency analysis.

READ-ONLY use of ../../qubo.py, ../../squeeze/mm.py, ../../squeeze/mmqb.py.

The QUBO built by mm.build_modmul is, structurally:
    variables x in {0,1}^n  (input bits, word bits, AND vars, carry/adder bits)
    for every add_square(lin, const):   lin . x + const == 0     (linear equality)
    for every AND cache entry z=(i,j):   z == x_i AND x_j          (nonconvex)
E == 0  iff  all equalities hold and all AND vars are consistent.

LP relaxation P:
    box            0 <= x <= 1
    equalities     lin . x + const == 0          (exact, kept as equalities)
    McCormick(AND) z<=x_i, z<=x_j, z>=x_i+x_j-1, z>=0   (convex hull of z=x_i x_j)
P contains the integer feasible set, so any coordinate constant over P is
constant over every integer solution (sound removal).
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'squeeze'))

from mmqb import MMQB, Wrd     # noqa: E402
import mm                      # noqa: E402


def build_modmul_instance(p, mult='schoolbook', red='quotient', leaf=32, mode='wallace',
                          square=False):
    """Build A*B == C (mod p) with A,B as free input words. Returns dict with the
    compiler Q, the input Wrds, and the output word."""
    Q = MMQB(chunk=16, mode=mode)
    s = p.bit_length()

    def mkinput(name, nb):
        bits = [Q.new(f"{name}[{t}]", 'input') for t in range(nb)]
        Q.n_word += nb
        Q.trace.append(('word', name, bits, (lambda wv, name=name: wv[name])))
        return Wrd(bits, (lambda wv, name=name: wv[name]), name)

    A = mkinput("A", s)
    B = A if square else mkinput("B", s)
    Cw = Q.mkword("C", s, lambda wv: (wv["A"] * wv["B"]) % p)
    mm.build_modmul(Q, p, A, B, Cw, mult=mult, red=red, leaf=leaf, tag='mm')
    Q.finalize()
    return dict(Q=Q, p=p, s=s, A=A, B=B, C=Cw, square=square)


class LP:
    """LP relaxation of an instance's constraint system."""
    def __init__(self, inst):
        Q = inst['Q']
        self.Q = Q
        self.n = Q.n
        self.inst = inst
        # equalities from squares:  lin.x + const == 0
        self.eq = [(dict(lin), const) for lin, const in Q.squares]
        # AND relations
        self.ands = [(z, i, j) for (i, j), z in Q.andcache.items()]
        self.input_bits = [v for v in range(Q.n) if Q.kind[v] == 'input']
        self._build_matrices()

    def _build_matrices(self):
        n = self.n
        # equality rows
        A_eq_rows, b_eq = [], []
        for lin, const in self.eq:
            row = np.zeros(n)
            for v, w in lin.items():
                row[v] += w
            A_eq_rows.append(row)
            b_eq.append(-const)         # lin.x = -const
        # McCormick inequalities  A_ub x <= b_ub
        A_ub_rows, b_ub = [], []
        for z, i, j in self.ands:
            # z - x_i <= 0
            r = np.zeros(n); r[z] = 1; r[i] = -1; A_ub_rows.append(r); b_ub.append(0.0)
            r = np.zeros(n); r[z] = 1; r[j] = -1; A_ub_rows.append(r); b_ub.append(0.0)
            # -z + x_i + x_j <= 1
            r = np.zeros(n); r[z] = -1; r[i] = 1; r[j] = 1; A_ub_rows.append(r); b_ub.append(1.0)
            # z >= 0 handled by box
        self.A_eq = np.array(A_eq_rows) if A_eq_rows else np.zeros((0, n))
        self.b_eq = np.array(b_eq) if b_eq else np.zeros(0)
        self.A_ub = np.array(A_ub_rows) if A_ub_rows else np.zeros((0, n))
        self.b_ub = np.array(b_ub) if b_ub else np.zeros(0)

    def bounds(self, fixed=None):
        fixed = fixed or {}
        bnds = []
        for v in range(self.n):
            if v in fixed:
                bnds.append((float(fixed[v]), float(fixed[v])))
            else:
                bnds.append((0.0, 1.0))
        return bnds

    def solve(self, c, fixed=None):
        from scipy.optimize import linprog
        res = linprog(c, A_ub=self.A_ub if self.A_ub.shape[0] else None,
                      b_ub=self.b_ub if self.b_ub.shape[0] else None,
                      A_eq=self.A_eq if self.A_eq.shape[0] else None,
                      b_eq=self.b_eq if self.b_eq.shape[0] else None,
                      bounds=self.bounds(fixed), method='highs')
        return res

    def persistency(self, fixed=None, tol=1e-6, n_probe=6, seed=0, verbose=False):
        """Return dict v-> fixed_value for every coordinate constant over P.
        Filter with random-objective solves, confirm survivors with min/max."""
        rng = np.random.default_rng(seed)
        n = self.n
        # baseline feasibility
        base = self.solve(np.zeros(n), fixed)
        if not base.success:
            raise RuntimeError("LP infeasible: " + base.message)
        val = {v: base.x[v] for v in range(n)}
        lo = base.x.copy(); hi = base.x.copy()
        for _ in range(n_probe):
            c = rng.standard_normal(n)
            r = self.solve(c, fixed)
            if r.success:
                lo = np.minimum(lo, r.x); hi = np.maximum(hi, r.x)
            r = self.solve(-c, fixed)
            if r.success:
                lo = np.minimum(lo, r.x); hi = np.maximum(hi, r.x)
        candidates = [v for v in range(n) if hi[v] - lo[v] <= tol]
        # confirm each candidate with true min/max LP
        fixed_vars = {}
        for v in candidates:
            c = np.zeros(n); c[v] = 1.0
            rmin = self.solve(c, fixed)
            rmax = self.solve(-c, fixed)
            if rmin.success and rmax.success:
                vmin, vmax = rmin.x[v], rmax.x[v]
                if vmax - vmin <= tol:
                    fixed_vars[v] = 0.5 * (vmin + vmax)
        return fixed_vars, lo, hi


def witness_full(inst, a, b):
    """integer solution for given inputs; returns x array."""
    Q = inst['Q']
    wv0 = {"A": a}
    if not inst['square']:
        wv0["B"] = b
    x, wv = Q.witness({}, wv0)
    return np.array(x), wv


if __name__ == '__main__':
    inst = build_modmul_instance(13)
    lp = LP(inst)
    print("n", lp.n, "eq", len(lp.eq), "and", len(lp.ands), "inputs", len(lp.input_bits))
    fx, lo, hi = lp.persistency()
    print("LP-fixed (inputs free):", len(fx), "/", lp.n)
