"""Full 2-bit scan over all 1156 boolean free inputs. Sharded.
usage: python3 sa/pair.py <shard> <nshards> <align:0|1>
"""
import sys, time, pickle, itertools, json, os
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/sa')
import lib

shard = int(sys.argv[1]); nsh = int(sys.argv[2]); al = bool(int(sys.argv[3]))
lib.init_base()
bits = lib.bfree
n = len(bits)
THRESH = 11

best = []
cnt = 0
t0 = time.time()
hits = []
for i in range(n):
    if i % nsh != shard:
        continue
    for j in range(i + 1, n):
        s, f, v, nz = lib.score([bits[i], bits[j]], alignment=al)
        cnt += 1
        if s <= THRESH:
            hits.append((s, bits[i], bits[j]))
        best.append((s, bits[i], bits[j]))
    if len(best) > 200000:
        best.sort(); best = best[:2000]
best.sort(); best = best[:2000]
tag = f'{"A" if al else "N"}{shard}'
pickle.dump({'best': best, 'hits': hits, 'count': cnt},
            open(f'sa/pair_{tag}.pkl', 'wb'))
lo = min((b[0] for b in best), default=None)
print(f'shard {tag}: {cnt} pairs in {time.time()-t0:.0f}s  min={lo}  n_le_{THRESH}={len(hits)}')
