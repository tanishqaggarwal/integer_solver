#!/usr/bin/env python3
"""optimality.py -- fit the comb cost law and locate its provable floor.

Cost of a comb encoding with A group additions and one-hot tables of size D = 2^w:

    cost(w) = ceil(256/w) * ( 3*M + mu(D) )

  M     = cost of one 256x256 modular multiplication in the QUBO
  mu(D) = cost of one one-hot table look-up of a 256-bit constant from D entries

Both are recovered here by least squares from the MEASURED per-window costs.
"""
import json, math

W = json.load(open('window256_neq.json'))


def fit(mode):
    pts = [(int(k.rsplit('_w', 1)[1]), v['vars']) for k, v in W.items()
           if k.startswith(mode + '_w')]
    pts.sort()
    # per_window(w) = a + b * 2^w        (a = 3M + fixed linear glue)
    xs = [1 << w for w, _ in pts]; ys = [v for _, v in pts]
    nn = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    b = (nn * sxy - sx * sy) / (nn * sxx - sx * sx)
    a = (sy - b * sx) / nn
    return a, b, pts


print("=" * 84)
print("COMB COST LAW, FITTED TO THE MEASURED 256-BIT WINDOWS")
print("=" * 84)
for mode in ('binary', 'wallace'):
    a, b, pts = fit(mode)
    print(f"\n{mode}:  per_window(w) = {a:,.0f} + {b:.1f} * 2^w")
    print(f"  => one modular multiplication M ~ {a/3:,.0f} qubits")
    print(f"  => one table entry            ~ {b:.1f} qubits")
    print(f"  {'w':>3} {'measured':>12} {'law':>12} {'err':>7} {'TOTAL':>14}")
    for w, v in pts:
        pred = a + b * (1 << w)
        M = math.ceil(256 / w)
        print(f"  {w:3d} {v:12,d} {pred:12,.0f} {100*(pred-v)/v:6.1f}% {M*v:14,d}")
    best = min(pts, key=lambda t: math.ceil(256 / t[0]) * t[1])
    print(f"  optimum w={best[0]}: {math.ceil(256/best[0])*best[1]:,d} qubits")

print()
print("=" * 84)
print("THE FLOOR")
print("=" * 84)
a, b, _ = fit('binary')
M = a / 3
print(f"""
Coverage bound.  A comb with A additions and tables of size D offers D^(A+1)
distinct digit tuples; to address every k in [0,2^256) we need

        D^(A+1) >= 2^256      i.e.      A >= 256/log2(D) - 1 = 256/w - 1.

So the additions cannot be traded away by widening the window faster than
logarithmically, while the look-up cost grows as 2^w.  That is exactly the
measured trade-off above, and it is why the optimum sits at w = 8-10 rather
than at either extreme.

Additions cost 3 multiplications each and no fewer: the affine group law needs
lam*(x2-x1) = y2-y1, lam^2 = x3+x1+x2, lam*(x1-x3) = y3+y1.  Three products,
one of them a squaring.  Projective coordinates trade the inverse for MORE
products, not fewer.

Therefore, within the comb family,

        cost(w) >= (256/w) * 3M     with M = {M:,.0f} measured.

For that to fall below a 4,400-qubit annealer we would need

        256/w <= 4400 / (3*{M:,.0f})  =>  w >= {256/(4400/(3*M)):,.0f}

i.e. a look-up table of 2^{256/(4400/(3*M)):,.0f} precomputed points.  Not a
close call.

The floor is therefore ONE modular multiplication, {M:,.0f} qubits, about {M/4400:.0f}x a
real annealer -- and no encoding of this decision problem can contain fewer
than one, because deciding it requires comparing x(kG) with x(T) in F_p.
""")
