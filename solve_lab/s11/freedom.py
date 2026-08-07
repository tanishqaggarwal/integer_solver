"""For each control: perturb it by 1 and ask whether EVERY check it sits in can be re-closed
   by a handle other than the control itself. If yes -> genuinely free mod p."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
BITS = (542, 47, 438, 91)
CTRL = [5096, 21589, 14515, 19750, 33708, 31339, 16441, 22917, 13222, 14681, 28486, 38667]

base = [0] * L.NVARS
for b in BITS:
    base[b] = 1
fw.forward(base)
BASEBAD = set(fw.bad_checks(base))
print("base bad:", len(BASEBAD), sorted(BASEBAD))

for c in CTRL:
    v = [x for x in base]
    v[c] = v[c] + 1
    fw.forward(v)
    newbad = [a for a in fw.bad_checks(v) if a not in BASEBAD]
    verdict = []
    for a in newbad:
        try:
            hs, bs = deep.handles(v, a, locked={c} | set(BITS))
        except Exception:
            hs, bs = [], 1
        ok = any(d and bs % d == 0 for _, d in hs)
        verdict.append((a, 'closable' if ok else 'STUCK'))
    free = all(t == 'closable' for _, t in verdict)
    print(f"x{c}: perturbation breaks {len(newbad)} checks -> {verdict}   => {'FREE' if free else 'PINNED'}", flush=True)
