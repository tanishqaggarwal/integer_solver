"""The six coupled mod-p residuals in channel U=0,V=1, and an exhaustive control scan."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ['x3719', 'x25118', 'a14445', 'a34580', 'a33796', 'a27139']
TV = (3719, 25118, 33129, 3757, 33708, 10170, 31339, 6858, 37088, 13585, 12000)

CS = set()
st = list(TV)
while st:
    u = st.pop()
    if u in CS:
        continue
    CS.add(u)
    dd = L.definer.get(u)
    if dd is None:
        continue
    for w in L.avars[dd]:
        if w != u:
            st.append(w)
ORDER = [c for c in fw.ORDER if any(u in CS for u in c)]
print(f"[sys6] cone {len(CS)} vars, {len(ORDER)} comps of {len(fw.ORDER)}", file=sys.stderr)

BASE = [0] * L.NVARS
d = json.load(open(os.path.join(HERE, 'data', 'three.json')))
for k, x in d.items():
    BASE[int(k)] = int(x)


def lfwd(v):
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
    return lfwd(v)


def six(v):
    return [v[3719] % P, v[25118] % P, (v[33129] - v[3757]) % P, (v[33708] - v[10170]) % P,
            (v[31339] - v[6858]) % P, (v[37088] - v[13585]) % P]


if __name__ == '__main__':
    import random
    rnd = random.Random(5)
    FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None]
    # generic point: randomise the cone's free inputs so no derivative vanishes by accident
    cone_free = [u for u in sorted(CS) if L.definer.get(u) is None]
    th = {u: BASE[u] + rnd.randrange(1, 1 << 40) for u in cone_free}
    base = six(ev(th))
    print("residual names:", NAMES)
    t0 = time.time()
    rows = {i: [] for i in range(6)}
    for n, c in enumerate(FREE):
        t2 = dict(th)
        t2[c] = th.get(c, BASE[c]) + 1
        r1 = six(ev(t2))
        for i in range(6):
            if (r1[i] - base[i]) % P:
                rows[i].append(c)
        if n % 2000 == 0:
            print(f"  {n}/{len(FREE)} ({time.time()-t0:.0f}s)", flush=True)
    print(f"scan done ({time.time()-t0:.0f}s)")
    for i in range(6):
        print(f"  {NAMES[i]:8s}: {len(rows[i])} -> {rows[i][:16]}")
    json.dump({NAMES[i]: rows[i] for i in range(6)}, open(os.path.join(HERE, 'data', 'sys6.json'), 'w'))
