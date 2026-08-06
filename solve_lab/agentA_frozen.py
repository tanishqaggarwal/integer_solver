#!/usr/bin/env python3
"""Determine code dimension: which bits are 'frozen' (sole breaker of some check).
Also compute GF(2) null space of the single-bit incidence matrix."""
import json
from collections import defaultdict

cd = json.load(open('agentA_code.json'))
allbits = cd['allbits']
perbit = {int(k): set(v) for k, v in cd['perbit'].items()}

# For each check, which bits break it (single-flip)
check_bits = defaultdict(set)
for b in allbits:
    for c in perbit[b]:
        check_bits[c].add(b)

sole = {}   # check -> the unique bit that breaks it
for c, bs in check_bits.items():
    if len(bs) == 1:
        sole[c] = next(iter(bs))
frozen = set(sole.values())
print(f"total bits: {len(allbits)}")
print(f"checks broken by exactly 1 bit: {sum(1 for c in check_bits if len(check_bits[c])==1)}")
print(f"FROZEN bits (sole breaker of >=1 check): {len(frozen)}")
free_bits = [b for b in allbits if b not in frozen]
print(f"potentially-free bits: {len(free_bits)} -> {free_bits}")

# GF(2) incidence matrix null space (checks x bits)
checks = sorted(check_bits)
cidx = {c: i for i, c in enumerate(checks)}
# rows = checks as bitmask over bits
bidx = {b: i for i, b in enumerate(allbits)}
rows = []
for c in checks:
    m = 0
    for b in check_bits[c]:
        m |= (1 << bidx[b])
    rows.append(m)

# Gaussian elimination over GF(2) to find rank and nullity (num bits = 78 cols)
ncol = len(allbits)
basis = []   # pivot rows
pivcol = []
mat = rows[:]
r = 0
# reduce
pivots = {}
for m in mat:
    cur = m
    for pc, pm in pivots.items():
        if (cur >> pc) & 1:
            cur ^= pm
    if cur:
        # find lowest set bit as pivot
        lb = (cur & -cur).bit_length() - 1
        pivots[lb] = cur
rank = len(pivots)
nullity = ncol - rank
print(f"\nGF(2) incidence: {len(checks)} checks x {ncol} bits, rank={rank}, nullity={nullity}")
print(f"=> GF(2) code dimension (null space of single-bit incidence) = {nullity}")

json.dump({'frozen': sorted(frozen), 'free_bits': free_bits,
           'gf2_rank': rank, 'gf2_nullity': nullity}, open('agentA_frozen.json', 'w'))
print("saved agentA_frozen.json")
