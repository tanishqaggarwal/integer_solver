import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast, leaves
P = L.P

# control cone: free vars feeding the six targets
CONE = set()
for t in (3719, 25118, 25614, 34220, 12186, 1308, 24908, 19083):
    fr, _ = leaves.cone_free(t)
    CONE |= fr
# exclude boolean message bits (cones of the four activation trees)
BOOL = set()
for t in (8599, 21839, 7304, 25956):
    fr, _ = leaves.cone_free(t)
    BOOL |= fr
C = sorted(CONE - BOOL)
print(f"control cone: {len(CONE)} free, minus {len(CONE & BOOL)} boolean -> {len(C)} controls")

if __name__ == '__main__':
    rnd = random.Random(11)
    th = {c: rnd.randrange(1, 1 << 40) for c in C}
    base = fast.targets(fast.light(th))
    t0 = time.time()
    rows = {i: [] for i in range(6)}
    for c in C:
        t2 = dict(th)
        t2[c] = th[c] + 1
        r1 = fast.targets(fast.light(t2))
        for i in range(6):
            d = (r1[i] - base[i]) % P
            if d:
                rows[i].append(c)
    print(f"({time.time()-t0:.0f}s)  at a GENERIC point:")
    for i in range(6):
        print(f"  {fast.NAMES[i]:8s}: {len(rows[i])} controls -> {rows[i][:14]}")
    json.dump(C, open('cone_controls.json', 'w'))
