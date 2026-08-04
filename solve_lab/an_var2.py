#!/usr/bin/env python3
import re,json
VAR_RE=re.compile(r'x_(\d+)')
lines=open('../EQUATIONS.txt').read().split('\n')
# free inputs: read from forward_construct if it dumps, else infer.
# Find all equations containing x_4432 and x_7068, classify.
for tgt in [4432,7068]:
    eqs=[i for i,L in enumerate(lines) if L.strip() and f"x_{tgt}" in set(f"x_{m}" for m in VAR_RE.findall(L))]
    print(f"=== x_{tgt}: appears in {len(eqs)} eqs: {eqs[:30]}")
# Check forward_construct free-input detection
print("\n=== forward_construct.py head (free input logic) ===")
