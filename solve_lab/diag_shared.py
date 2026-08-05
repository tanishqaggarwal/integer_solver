#!/usr/bin/env python3
"""Look at the linear equations shared by the remainder free inputs. Determine whether they pin
the load residue mod p, or leave freedom. Show the structure of eq 1679, 4028, 4279, 4329."""
import json, re, ast, sys
from propagate import NVARS
p=2**256-2**32-977
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
for i in [1679,4028,4279,4329,5564]:
    L=lines[i].rsplit('=',1)[0]
    vs=sorted(set(int(m) for m in VAR.findall(L)))
    print(f"\n=== eq {i} ({len(L)} chars), vars={vs}")
    print("   ", L[:400])
