"""Cone-restricted fast evaluation of the six structural targets."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
BITS = (542, 47, 438, 91)
TVARS = (3719, 25118, 25614, 34220, 12186, 1308, 24908, 19083)

# variable cone of the targets
CS = set()
st = list(TVARS)
while st:
    u = st.pop()
    if u in CS:
        continue
    CS.add(u)
    d = L.definer.get(u)
    if d is None:
        continue
    for w in L.avars[d]:
        if w != u:
            st.append(w)
ORDER = [c for c in fw.ORDER if any(u in CS for u in c)]
NEED = [u for u in CS if L.definer.get(u) is not None]
print(f"[quick] cone {len(CS)} vars, {len(ORDER)} scc comps (of {len(fw.ORDER)})", file=sys.stderr)

BASE = [0] * L.NVARS
for b in BITS:
    BASE[b] = 1


def lforward(v):
    for comp in ORDER:
        if len(comp) == 1:
            u = comp[0]
            x = fw.solve_lin(L.definer[u], u, v)
            if x is not None:
                v[u] = x
        else:
            for _ in range(40):
                ch = False
                for u in comp:
                    x = fw.solve_lin(L.definer[u], u, v)
                    if x is not None and x != v[u]:
                        v[u] = x
                        ch = True
                if not ch:
                    break
    return v


def ev(th):
    v = BASE[:]
    for k, x in th.items():
        v[k] = x
    return lforward(v)


def six(v):
    return [v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P,
            (v[12186] - v[1308]) % P, (v[24908] - v[19083]) % P]


if __name__ == '__main__':
    import random
    t0 = time.time()
    rnd = random.Random(1)
    th = {c: rnd.randrange(1, 1 << 60) for c in [16441, 22917, 31339, 33708, 5096, 21589, 14515, 19750]}
    for _ in range(20):
        ev(th)
    print(f"20 evals in {time.time()-t0:.2f}s")
    # cross-check against full forward
    import fast
    v1 = ev(th)
    v2 = fast.light(th)
    print("matches full forward:", six(v1) == six(v2))
