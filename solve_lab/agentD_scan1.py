#!/usr/bin/env python3
"""Scan 1: the 4 quadrants with default constant routing, plus timing.
Quadrant set by activators: (1,1)=a7+a34, (1,0)=a7 only, (0,1)=a34 only, (0,0)=none."""
import time, json
import agentD_harness as H
C1,C2=H.C1,H.C2
a7=24601; a34=2081
CONST={30213:C2, 22162:C1, 24468:C1, 18956:C2}
def show(name, r):
    print(f"{name:10s} sat={r['satisfied']} nfail={r['nfail']} x15298={r['x_15298']} "
          f"x7715={r['x_7715']} x34554={r['x_34554']} S0={r['S_is0']} T0={r['T_is0']} "
          f"core_fail={r['core_fail']} noncore={r['noncore_fail']} det={r['determined']}", flush=True)

t0=time.time()
r11=H.run_config({a7:1, a34:1, **CONST}); show("(1,1)", r11)
print(f"  [per-config {time.time()-t0:.1f}s]", flush=True)
r10=H.run_config({a7:1, **CONST}); show("(1,0)", r10)
r01=H.run_config({a34:1, **CONST}); show("(0,1)", r01)
r00=H.run_config({**CONST}); show("(0,0)+C", r00)
r00b=H.run_config({}); show("(0,0)noC", r00b)
print("F(1,0):", r10['F'][:30], flush=True)
print("F(0,1):", r01['F'][:30], flush=True)
print("F(0,0)+C:", r00['F'][:30], flush=True)
print("F(0,0)noC:", r00b['F'][:40], flush=True)
