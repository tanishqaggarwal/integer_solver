#!/usr/bin/env python3
"""Baby-step/giant-step for a SMALL scalar: is k (or N-k) below 2^44?
Also covers k = c*2^e for small c (the ladder is closed under doubling)."""
import json, time, sys
import numpy as np
from model import P, TARGET, to_short
from group import add, neg, mul
from fastgrp import batch_add

lad = json.load(open('ladder.json'))['ladder']
G = (int(lad[0][1]), int(lad[0][2]))
T = to_short(TARGET)
M = 1 << 22                      # baby steps -> covers k < 2^44
t0 = time.time()
keys = np.empty(M, dtype=np.uint64)
Q = None
pts = []
cur = None
for j in range(M):
    cur = G if cur is None else add(cur, G)
    keys[j] = cur[0] & ((1 << 63) - 1)
    if j % 500000 == 0: print(' baby', j, '%.0fs' % (time.time() - t0), flush=True)
ordr = np.argsort(keys); skeys = keys[ordr]
print('baby table built %.0fs' % (time.time() - t0), flush=True)
S = mul(M, G); nS = neg(S)
cur = T
found = None
for i in range(M):
    kk = np.uint64(cur[0] & ((1 << 63) - 1))
    idx = np.searchsorted(skeys, kk)
    while idx < M and skeys[idx] == kk:
        j = int(ordr[idx]) + 1
        k = i * M + j
        if mul(k, G) == T: found = k; break
        idx += 1
    if found: break
    cur = add(cur, nS)
    if i % 200000 == 0: print(' giant', i, '%.0fs' % (time.time() - t0), flush=True)
print('BSGS result:', found, '%.0fs' % (time.time() - t0), flush=True)
json.dump({'bound': 2 ** 44, 'k': str(found) if found else None}, open('runs/bsgs.json', 'w'))
