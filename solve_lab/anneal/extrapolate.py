#!/usr/bin/env python3
"""extrapolate.py -- from the measured cost model to the real 256-bit instance."""
import json, math

S = 256                     # field size of the real instance
MBITS = 256                 # scalar size

# measured:  vars per comb window  =  A*s^2 + (alpha*2^w + beta)*s
MODEL = {   # A,     alpha,  beta      (least-squares over s=16..64, w=1,4,8)
    'wallace': (13.70, 2.044, 75.1),
    'binary':  (5.65,  0.375, 45.25),
}
COUP = {    # couplers per window (same least-squares procedure)
    'wallace': (60.82, 9.716, 265.9),
    'binary':  (108.71, 20.311, 1097.8),
}


def per_window(mode, w, s=S, tbl=MODEL):
    A, al, be = tbl[mode]
    return A * s * s + (al * (1 << w) + be) * s


def total(mode, w, mbits=MBITS, s=S):
    return math.ceil(mbits / w) * per_window(mode, w, s)


def scan(mode):
    best = min(range(1, 17), key=lambda w: total(mode, w))
    return best, total(mode, best)


if __name__ == '__main__':
    print("=" * 78)
    print("MONOLITHIC ENCODING OF THE FULL INSTANCE  (256-bit field, 256-bit scalar)")
    print("=" * 78)
    print(f"{'mode':>8} {'w':>3} {'windows':>8} {'logical qubits':>16} {'couplers':>14}")
    for mode in ('binary', 'wallace'):
        for w in (1, 4, 6, 8, 10, 12):
            M = math.ceil(MBITS / w)
            print(f"{mode:>8} {w:3d} {M:8d} {M*per_window(mode,w):16,.0f} "
                  f"{M*per_window(mode,w,tbl=COUP):14,.0f}")
        w, t = scan(mode)
        print(f"  --> optimum w={w}: {t:,.0f} logical qubits\n")

    print("=" * 78)
    print("SMALLEST INDIVISIBLE PIECES")
    print("=" * 78)
    for mode in ('binary', 'wallace'):
        A = MODEL[mode][0]
        print(f"{mode:>8}: one 256x256 modular multiplication  ~ {A/4*S*S:12,.0f} qubits")
        print(f"{mode:>8}: one comb window (4 modmuls + MUX)   ~ {per_window(mode,8):12,.0f} qubits")
    print()

    print("=" * 78)
    print("WHAT FITS ON REAL HARDWARE")
    print("=" * 78)
    for name, N, deg in (("D-Wave Advantage  (Pegasus)", 5760, 15),
                         ("D-Wave Advantage2 (Zephyr)", 4400, 20)):
        print(f"\n{name}: {N} qubits, degree {deg}")
        for mode in ('binary', 'wallace'):
            A = MODEL[mode][0]
            s_mul = int((N / (A / 4)) ** 0.5)
            # largest complete instance: total(s, m=s, w) <= N
            best = 0, 0
            for w in range(1, 9):
                s = 4
                while total(mode, w, mbits=s, s=s) <= N: s += 1
                if s - 1 > best[0]: best = (s - 1, w)
            print(f"  {mode:>8}: largest single modmul  = {s_mul:3d}-bit field")
            print(f"  {mode:>8}: largest COMPLETE DLP   = {best[0]:3d}-bit curve (w={best[1]}), "
                  f"vs the {S}-bit instance")
        print(f"  qubits short of the monolithic encoding: "
              f"{scan('binary')[1]/N:,.0f}x (before minor-embedding overhead)")

    print()
    print("=" * 78)
    print("HOW MANY EXECUTIONS BUY WHAT  (split k = k_hi*2^mu + k_lo, k_hi enumerated")
    print("classically, only the mu low bits annealed; w = 8)")
    print("=" * 78)
    print(f"{'mu':>5} {'windows':>8} {'qubits/run':>14} {'runs':>12} {'verdict':>34}")
    for mu in (8, 16, 32, 64, 128, 192, 256):
        M = math.ceil(mu / 8)
        q = M * per_window('binary', 8)
        runs = 2 ** (MBITS - mu)
        v = ("outer loop is the whole search" if mu < 128 else
             "no better than Pollard rho (2^128)" if mu < 256 else "single shot")
        print(f"{mu:5d} {M:8d} {q:14,.0f} {('2^%d' % (MBITS-mu)):>12} {v:>34}")
