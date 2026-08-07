import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
P = L.P

FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None]
th = {}
base = fast.targets(fast.light(th))
t0 = time.time()
rows = {i: [] for i in range(6)}
for n, c in enumerate(FREE):
    t2 = {c: 1}
    r1 = fast.targets(fast.light(t2))
    for i in range(6):
        d = (r1[i] - base[i]) % P
        if d:
            rows[i].append((c, d))
    if n % 1000 == 0:
        print(f"  {n}/{len(FREE)} ({time.time()-t0:.0f}s)", flush=True)
print(f"scan done ({time.time()-t0:.0f}s)")
for i in range(6):
    print(f"  {fast.NAMES[i]:8s}: {len(rows[i])} controls -> {[c for c,_ in rows[i][:14]]}")
json.dump({fast.NAMES[i]: [[c, str(d)] for c, d in rows[i]] for i in range(6)},
          open('scanall.json', 'w'))
