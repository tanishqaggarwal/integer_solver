#!/usr/bin/env python3
"""Concatenate per-thread shard parts, sort each shard in place, delete parts.
Because shard id = key>>61 is order-consistent with key order, per-shard sorting
gives an exact binary-searchable table (searched per shard by aa_signed)."""
import numpy as np, os, sys, glob, time
pref = sys.argv[1]
t0 = time.time(); tot = 0
for sh in range(8):
    parts = sorted(glob.glob('%s.%d.*' % (pref, sh)))
    if not parts:
        print('shard %d: no parts' % sh); continue
    arrs = [np.fromfile(f, dtype=np.uint64) for f in parts]
    a = np.concatenate(arrs) if len(arrs) > 1 else arrs[0]
    del arrs
    assert len(a) == 0 or ((a >> np.uint64(61)) == np.uint64(sh)).all(), 'shard %d mis-bucketed' % sh
    a.sort(kind='stable')
    a.tofile('%s.%d' % (pref, sh))
    tot += len(a)
    print('shard %d: %d keys  %.1fs' % (sh, len(a), time.time() - t0)); sys.stdout.flush()
    del a
    for f in parts: os.remove(f)
print('TOTAL %d keys  %.1fs' % (tot, time.time() - t0))
