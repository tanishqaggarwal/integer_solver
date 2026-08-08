#!/usr/bin/env python3
"""synth_circuit.py -- a SAME-SHAPE synthetic instance of the full problem.

We compile a small planted-key ECDLP into an EQUATIONS.txt-style polynomial system
over integer unknowns x_i, using EXACTLY the gadget algebra that CIRCUIT_STRUCTURE.md
reverse-engineered from the real instance:

  * 256-selector-bit comb over a doubling chain  L_i = 2^i G   (here m bits)
  * an OFFSET running point  R_0 = Q0  (signed/offset digit bookkeeping) so the
    accumulator is never the identity O and never hits a doubling/degenerate add
  * per bit i a GATED, fraction-free COMBINE gadget computing R_i + L_i:
        A = (w5 + w1 + w3 + a2)*(w3 - w1)^2 - (w4 - w2)^2     ==  p * sA
        B = (w6 + w2)*(w3 - w1) - (w4 - w2)*(w1 - w5)         ==  p * sB
    with (w1,w2)=R_i, (w3,w4)=L_i, (w5,w6)=addpt_i, gated by the selector b_i
    through  b_i*(c1*A + c2*B) - p*free = 0  (the p-slack encoding of "= 0 mod p")
  * the d != 0 NON-DEGENERACY gadget that closes the x1 == x2 division loophole
    (CIRCUIT_STRUCTURE.md soundness audit item 1): a STRONG inverse wire
        b_i*(d_i * dinv_i - 1) - p * sD = 0 ,   d_i = w3 - w1
    so an active gadget forces d_i to be a unit mod p, i.e. w3 != w1 mod p.
  * a 4-way / here 2-way MUX routing R_{i+1} = b_i ? (R_i + L_i) : R_i
  * bit-gated constant LOADS of the leaf coordinates  b*(v-C)=0, (1-b)*v=0
  * boolean atoms  b*(b-1)=0  on every selector
  * a hardwired  a2  pin  (a2 - x_a2 = 0, a2 = 0 here)  used by every gadget
  * a 220-style COPY class linking several wires  (x_a - x_b = 0)
  * the final register-vs-target comparison
        Rroot.x - (Q0+T).x = 0 ,   Rroot.y - (Q0+T).y = 0
    plus OR(all selectors) = 1

The equations are FIXED.  The 256 (here m) selector bits are the ONLY free inputs.
reconstruct(bits) recomputes EVERY x_i deterministically by honest point arithmetic;
the resulting assignment satisfies every equation exactly in Z  <=>  sum b_i 2^i G = T.

Everything is exact big-integer arithmetic; nothing is mod-p-only (the p-slack wires
carry the exact quotients), so solve_lab/checker.py validates it verbatim.
"""
import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))  # solve_lab/..
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))         # anneal/
from synth.gen import make


class Builder:
    """Allocates x_i indices and accumulates (equation-string, planted-value) data."""
    def __init__(self):
        self.n = 0
        self.eqs = []          # list of equation LHS strings ( "... " meaning LHS = 0 )
        self.recipe = []       # list of ('kind', payload) to recompute values from bits
        self.planted = {}      # idx -> planted integer value (for self-test only)

    def var(self, val=0):
        i = self.n; self.n += 1
        self.planted[i] = int(val)
        return i

    def emit(self, lhs):
        self.eqs.append(lhs)


def _fit_offset(inst, seed_scan=64):
    """Pick an offset scalar q0 so the whole planted accumulation path is generic:
    every active add R + L_i has distinct x-coords and no intermediate point is O."""
    c, G, n, k, m = inst.curve, inst.G, inst.n, inst.k, inst.bits
    bits = [(k >> i) & 1 for i in range(m)]
    leaves = inst.pts(m)                    # L_i = 2^i G
    for q0 in range(3, 3 + seed_scan):
        R = c.mul(q0, G)
        ok = R is not None
        path = []
        for i in range(m):
            Li = leaves[i]
            if bits[i]:
                # active add R + Li : need R != O, Li != O, R.x != Li.x, R != -Li
                if R is None or Li is None or R[0] == Li[0]:
                    ok = False; break
                Rn = c.add(R, Li)
                if Rn is None: ok = False; break
                path.append(Rn); R = Rn
            else:
                path.append(R)
        if ok and R is not None and R == c.add(c.mul(q0, G), inst.T):
            return q0
    return None


def build(inst):
    """Compile inst into (equations, meta, reconstruct-recipe). Returns a dict."""
    c, G, n, k, m = inst.curve, inst.G, inst.n, inst.k, inst.bits
    p, Bc = c.p, c.B
    a2 = 0                                   # y^2 = x^3 + B  => a2 = 0 (short form)
    leaves = inst.pts(m)                     # constants L_i = 2^i G, i=0..m-1
    q0 = _fit_offset(inst)
    if q0 is None:
        raise RuntimeError("no generic offset found; reseed the instance")
    Q0 = c.mul(q0, G)
    Ttot = c.add(Q0, inst.T)                 # target the root must equal: Q0 + k G
    bplant = [(k >> i) & 1 for i in range(m)]

    B_ = Builder()
    # ---- selector bits: the ONLY free inputs ----
    bit_idx = [B_.var(bplant[i]) for i in range(m)]
    for i in range(m):
        # boolean atom  b*(b-1) = 0
        B_.emit(f"(x_{bit_idx[i]})*((x_{bit_idx[i]})-(1))")

    # ---- the a2 pin (hardwired constant gate, a2 = 0) ----
    a2_idx = B_.var(a2)
    B_.emit(f"(x_{a2_idx})-({a2})")

    # ---- per-bit leaf LOAD (bit-gated constant) ----
    lx_idx, ly_idx = [], []
    for i in range(m):
        Lx, Ly = leaves[i]
        xi = B_.var(Lx if bplant[i] else 0)
        yi = B_.var(Ly if bplant[i] else 0)
        b = bit_idx[i]
        # b*(v - C) = 0    and   (1-b)*v = 0
        B_.emit(f"(x_{b})*((x_{xi})-({Lx}))")
        B_.emit(f"(1-(x_{b}))*(x_{xi})")
        B_.emit(f"(x_{b})*((x_{yi})-({Ly}))")
        B_.emit(f"(1-(x_{b}))*(x_{yi})")
        lx_idx.append(xi); ly_idx.append(yi)

    # ---- offset root R_0 = Q0 (two pinned constant registers) ----
    rx = B_.var(Q0[0]); ry = B_.var(Q0[1])
    B_.emit(f"(x_{rx})-({Q0[0]})")
    B_.emit(f"(x_{ry})-({Q0[1]})")

    # random small combining coefficients for the gadget residuals (like the real 3 combos)
    import random
    rnd = random.Random(0xC0FFEE ^ m)

    R = Q0                                    # honest running point (for planting)
    stage_recipe = []                         # to recompute at reconstruction time
    for i in range(m):
        b = bit_idx[i]
        w1, w2 = rx, ry                       # R_i registers
        w3, w4 = lx_idx[i], ly_idx[i]         # L_i registers (gated-loaded)
        Rx_val, Ry_val = B_.planted[w1], B_.planted[w2]
        Lx, Ly = leaves[i]

        if bplant[i]:
            addpt = c.add(R, leaves[i])
            ax, ay = addpt
            d_val = (Lx - Rx_val) % p
            dinv_val = pow(d_val, -1, p)
        else:
            # gated OFF: gadget inactive. Put honest-but-inert values.
            ax, ay = 0, 0
            d_val = 0; dinv_val = 0

        w5 = B_.var(ax); w6 = B_.var(ay)      # addpt registers
        dinv = B_.var(dinv_val)
        # --- the two fraction-free residuals, p-slack encoded, gated by b ---
        # A = (w5+w1+w3+a2)*(w3-w1)^2 - (w4-w2)^2
        Aexpr = (f"((x_{w5})+(x_{w1})+(x_{w3})+(x_{a2_idx}))*"
                 f"((x_{w3})-(x_{w1}))*((x_{w3})-(x_{w1}))"
                 f"-((x_{w4})-(x_{w2}))*((x_{w4})-(x_{w2}))")
        # B = (w6+w2)*(w3-w1) - (w4-w2)*(w1-w5)
        Bexpr = (f"((x_{w6})+(x_{w2}))*((x_{w3})-(x_{w1}))"
                 f"-((x_{w4})-(x_{w2}))*((x_{w1})-(x_{w5}))")
        c1 = rnd.randrange(1, 40); c2 = rnd.randrange(1, 40)
        # value of A, B at planted assignment (exact integers)
        if bplant[i]:
            A_int = (ax + Rx_val + Lx + a2) * (Lx - Rx_val) ** 2 - (Ly - Ry_val) ** 2
            B_int = (ay + Ry_val) * (Lx - Rx_val) - (Ly - Ry_val) * (Rx_val - ax)
            combo = c1 * A_int + c2 * B_int           # divisible by p
            sA_val = combo // p
            assert combo % p == 0, "residual not divisible by p"
        else:
            sA_val = 0
        sA = B_.var(sA_val)
        # b*(c1*A + c2*B) - p*sA = 0
        B_.emit(f"(x_{b})*(({c1})*({Aexpr})+({c2})*({Bexpr}))-({p})*(x_{sA})")

        # --- non-degeneracy: b*(d*dinv - 1) - p*sD = 0 ,  d = w3 - w1 ---
        if bplant[i]:
            dprod = d_val_int = (Lx - Rx_val) * dinv_val - 1     # ≡ 0 mod p
            sD_val = ((Lx - Rx_val) * dinv_val - 1) // p
            assert ((Lx - Rx_val) * dinv_val - 1) % p == 0
        else:
            sD_val = 0
        sD = B_.var(sD_val)
        B_.emit(f"(x_{b})*(((x_{w3})-(x_{w1}))*(x_{dinv})-(1))-({p})*(x_{sD})")

        # --- 2-way MUX: R_{i+1} = b*(addpt) + (1-b)*R_i ---
        if bplant[i]:
            nxR, nyR = ax, ay
        else:
            nxR, nyR = Rx_val, Ry_val
        nrx = B_.var(nxR); nry = B_.var(nyR)
        B_.emit(f"(x_{nrx})-((x_{b})*(x_{w5})+(1-(x_{b}))*(x_{w1}))")
        B_.emit(f"(x_{nry})-((x_{b})*(x_{w6})+(1-(x_{b}))*(x_{w2}))")

        stage_recipe.append(dict(b=b, w3=w3, w4=w4, w5=w5, w6=w6, dinv=dinv,
                                 sA=sA, sD=sD, nrx=nrx, nry=nry, w1=w1, w2=w2,
                                 c1=c1, c2=c2, Lx=Lx, Ly=Ly))
        rx, ry = nrx, nry
        if bplant[i]:
            R = c.add(R, leaves[i])

    # ---- final comparison: root == Q0 + T ----
    B_.emit(f"(x_{rx})-({Ttot[0]})")
    B_.emit(f"(x_{ry})-({Ttot[1]})")

    # ---- OR(all selectors) = 1  (at least one bit set) : product of (1-b_i) = 0 ----
    prod = "*".join(f"(1-(x_{bit_idx[i]}))" for i in range(m))
    B_.emit(f"({prod})")

    # ---- a 220-style copy class: chain a few spare wires equal to a2 pin ----
    copy_prev = a2_idx
    copy_chain = []
    for _ in range(4):
        cvar = B_.var(B_.planted[a2_idx])
        B_.emit(f"(x_{cvar})-(x_{copy_prev})")
        copy_chain.append(cvar); copy_prev = cvar

    meta = dict(p=p, B=Bc, a2=a2, m=m, q0=q0, Q0=list(Q0), Ttot=list(Ttot),
                G=list(G), T=list(inst.T), k=k, n=n,
                bit_idx=bit_idx, a2_idx=a2_idx, lx_idx=lx_idx, ly_idx=ly_idx,
                rx0=None, ry0=None, root_x=rx, root_y=ry,
                leaves=[list(L) for L in leaves], nvars=B_.n,
                stage=stage_recipe, copy_chain=copy_chain, a2_pin=a2_idx)
    return dict(eqs=B_.eqs, meta=meta, planted=B_.planted)


def reconstruct(meta, bits, curve=None):
    """Deterministically recompute EVERY x_i from ONLY the selector bits.

    This is the Link-C reconstruction map: 256 (here m) bits -> full assignment.
    Returns a list v of length meta['nvars'].
    """
    from synth.gen import Curve
    p = meta['p']; a2 = meta['a2']; m = meta['m']
    c = curve or Curve(p, meta['B'])
    v = [0] * meta['nvars']
    # selectors
    for i, bi in enumerate(meta['bit_idx']):
        v[bi] = int(bits[i])
    # a2 pin
    v[meta['a2_idx']] = a2
    # leaf loads
    for i in range(m):
        Lx, Ly = meta['leaves'][i]
        v[meta['lx_idx'][i]] = Lx if bits[i] else 0
        v[meta['ly_idx'][i]] = Ly if bits[i] else 0
    # running point starts at Q0
    R = tuple(meta['Q0'])
    Rx, Ry = R
    for i, st in enumerate(meta['stage']):
        v[st['w1']]  # already R registers (rx/ry from previous stage or Q0 pins)
        Lx, Ly = st['Lx'], st['Ly']
        if bits[i]:
            add = c.add((Rx, Ry), (Lx, Ly))
            ax, ay = add
            d = (Lx - Rx) % p
            dinv = pow(d, -1, p)
            A_int = (ax + Rx + Lx + a2) * (Lx - Rx) ** 2 - (Ly - Ry) ** 2
            B_int = (ay + Ry) * (Lx - Rx) - (Ly - Ry) * (Rx - ax)
            sA = (st['c1'] * A_int + st['c2'] * B_int) // p
            sD = ((Lx - Rx) * dinv - 1) // p
            nx, ny = ax, ay
        else:
            ax = ay = 0; dinv = 0; sA = 0; sD = 0
            nx, ny = Rx, Ry
        v[st['w5']] = ax; v[st['w6']] = ay
        v[st['dinv']] = dinv; v[st['sA']] = sA; v[st['sD']] = sD
        v[st['nrx']] = nx; v[st['nry']] = ny
        Rx, Ry = nx, ny
    # copy class
    for cv in meta['copy_chain']:
        v[cv] = a2
    # the Q0 pin registers (rx,ry of stage 0's w1/w2) — set from Q0
    # (they are the first nrx/nry chain seed; find them: stage[0] w1,w2)
    if meta['stage']:
        v[meta['stage'][0]['w1']] = meta['Q0'][0]
        v[meta['stage'][0]['w2']] = meta['Q0'][1]
    return v


def write_equations(eqs, path):
    with open(path, 'w') as f:
        for lhs in eqs:
            f.write(lhs + " = 0\n")


if __name__ == '__main__':
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    inst = make(bits, seed=seed)
    art = build(inst)
    meta, eqs = art['meta'], art['eqs']
    print(f"synthetic same-shape instance: {bits}-bit key, p={inst.curve.p}")
    print(f"  planted k = {inst.k}")
    print(f"  variables x_0..x_{meta['nvars']-1}   ({meta['nvars']} vars)")
    print(f"  equations = {len(eqs)}")
    # self-test: planted assignment satisfies everything
    import re as _re
    vp = [art['planted'][i] for i in range(meta['nvars'])]
    codes = [compile(_re.sub(r'x_(\d+)', r'v[\1]', lhs), '<eq>', 'eval') for lhs in eqs]
    bad = [j for j, code in enumerate(codes) if eval(code, {'v': vp}) != 0]
    print(f"  planted assignment satisfies {len(eqs)-len(bad)}/{len(eqs)}"
          + (f"  FAIL@{bad[:5]}" if bad else "  (all)"))
    # reconstruct from bits and compare
    kb = [(inst.k >> i) & 1 for i in range(bits)]
    vr = reconstruct(meta, kb, inst.curve)
    print(f"  reconstruct(bits) == planted : {vr == vp}")
