#!/usr/bin/env python3
"""report.py -- headline resource table for the real instance, from measured window costs."""
import json, math

S = MB = 256
W = json.load(open('window256_neq.json'))
HW = [("D-Wave Advantage  (Pegasus, deg 15)", 5760),
      ("D-Wave Advantage2 (Zephyr,  deg 20)", 4400)]


def rows(mode):
    out = []
    for k, v in W.items():
        m, w = k.rsplit('_w', 1)
        if m != mode: continue
        w = int(w); M = math.ceil(MB / w)
        out.append((w, M, M * v['vars'], M * v['couplers'], v['jbits']))
    return sorted(out)


print("=" * 92)
print("MINIMAL MONOLITHIC QUBO FOR THE FULL INSTANCE")
print("  (256-bit field, 256-bit scalar; 3-multiplication affine step; measured, not modelled)")
print("=" * 92)
print(f"{'mode':>8} {'w':>3} {'windows':>8} {'qubits/window':>14} {'TOTAL qubits':>15} "
      f"{'TOTAL couplers':>16} {'|J| range':>10}")
best = {}
for mode in ('binary', 'wallace'):
    for w, M, tv, tc, jb in rows(mode):
        print(f"{mode:>8} {w:3d} {M:8d} {tv//M:14,d} {tv:15,d} {tc:16,d} {'2^%d'%jb:>10}")
        if mode not in best or tv < best[mode][0]: best[mode] = (tv, tc, w, jb)
    print()
print("OPTIMA")
for mode in ('binary', 'wallace'):
    tv, tc, w, jb = best[mode]
    print(f"  {mode:>8}: {tv:,d} qubits, {tc:,d} couplers at w={w}, coupler dynamic range 2^{jb}")
print()
print("=" * 92)
print("VERSUS REAL HARDWARE")
print("=" * 92)
for name, N in HW:
    tv, tc, w, jb = best['wallace']
    tvb = best['binary'][0]
    print(f"{name}: {N:,d} qubits")
    print(f"   shortfall, precision-safe encoding : {tv/N:12,.0f}x  ({tv:,d} / {N:,d})")
    print(f"   shortfall, qubit-minimal encoding  : {tvb/N:12,.0f}x  (but needs 2^21 coupler precision")
    print(f"                                          against ~4-5 bits available)")
