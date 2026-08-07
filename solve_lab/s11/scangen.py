import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
P = L.P
FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None]
BOOLSET = set()
import leaves
for t in (8599, 21839, 7304, 25956):
    fr, _ = leaves.cone_free(t)
    BOOLSET |= fr
rnd = random.Random(2026)
CONE = json.load(open('cone_controls.json'))
th = {c: rnd.randrange(1, 1 << 80) for c in CONE}
base = fast.targets(fast.light(th))
t0 = time.time()
rows = {i: [] for i in range(6)}
for n, c in enumerate(FREE):
    if c in BOOLSET:
        continue
    t2 = dict(th)
    t2[c] = th.get(c, 0) + 1
    r1 = fast.targets(fast.light(t2))
    for i in range(6):
        if (r1[i] - base[i]) % P:
            rows[i].append(c)
    if n % 1500 == 0:
        print(f"  {n}/{len(FREE)} ({time.time()-t0:.0f}s)", flush=True)
print(f"done ({time.time()-t0:.0f}s) -- GENERIC point, all free inputs:")
for i in range(6):
    print(f"  {fast.NAMES[i]:8s}: {len(rows[i])} -> {rows[i][:20]}")
json.dump({fast.NAMES[i]: rows[i] for i in range(6)}, open('scangen.json', 'w'))
