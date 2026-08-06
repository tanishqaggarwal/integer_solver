"""Hard-ordered closure: structural solve, then close each residual class with its own handle."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
BITS = (542, 47, 438, 91)
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
TH = {int(k): x for k, x in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','theta_solveB.json'))).items()}

v = [0] * L.NVARS
for b in BITS:
    v[b] = 1
for k, x in TH.items():
    v[k] = x
fw.forward(v)
print("six targets:", [v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P,
                       (v[12186]-v[1308]) % P, (v[24908]-v[19083]) % P])
v[14853] = v[12186]
v[16742] = v[24908]
v[30213] = G0
v[22820] = 0
v[7497] = (C0 + MM * G0) // P
v[22162] = -C0B
v[14393] = 0
v[11436] = 0
fw.forward(v)
LOCK = set(BITS) | set(TH) | {14853, 16742, 30213, 22820, 7497, 22162, 14393, 11436}
used = set()


def close_one(a, extra_lock=()):
    """close atom a with the cheapest available handle (shallow then deep)."""
    if fw.evalpoly(L.polys[a], v) == 0:
        return True
    lk = LOCK | used | set(extra_lock)
    cands = []
    for u in L.avars[a]:
        if L.definer.get(u) is None and u not in lk and not any(mm.count(u) > 1 for mm in L.polys[a]):
            cands.append((u, None))
    cands.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
    try:
        hs, base = deep.handles(v, a, locked=lk)
        cands += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
    except Exception:
        pass
    for t, d in cands:
        old = v[t]
        if d is None:
            x = fw.solve_lin(a, t, v)
            if x is None or x == old:
                continue
        else:
            bs = fw.evalpoly(L.polys[a], v)
            if not d or bs % d:
                continue
            x = old - bs // d
        v[t] = x
        fw.forward(v)
        if fw.evalpoly(L.polys[a], v) == 0:
            used.add(t)
            return True
        v[t] = old
        fw.forward(v)
    return False


STAGES = [
    ("group1 mirror", [26719, 26721, 26723]),
    ("group2 mirror", [26733, 28438, 32342]),
    ("linking pins", [21050, 34580, 33796, 15030, 9193, 14445, 35374, 26727, 31938, 31940]),
    ("load pins", [13438, 13440, 23824, 23826, 36040, 36042, 20103, 31986, 38567]),
    ("big 1-eq", [14312, 40047, 41774, 25676, 25793, 42245, 36185, 40812, 22932, 39719, 39760]),
]
t0 = time.time()
for name, atoms in STAGES:
    got = 0
    for a in atoms:
        if a >= L.NA:
            continue
        if close_one(a):
            got += 1
    bad = fw.bad_checks(v)
    print(f"  stage {name}: closed {got}/{len(atoms)} | bad={len(bad)} ({time.time()-t0:.0f}s) {bad[:12]}", flush=True)

# final sweep on whatever remains
for rnd in range(15):
    bad = fw.bad_checks(v)
    if not bad:
        break
    prog = False
    for a in sorted(bad, key=lambda a: len(L.atom2eq.get(a, {}))):
        if close_one(a):
            prog = True
    nb = fw.bad_checks(v)
    print(f"  sweep{rnd}: bad={len(nb)} ({time.time()-t0:.0f}s) {nb[:12]}", flush=True)
    if not prog or set(nb) == set(bad):
        break
fw.forward(v)
b = fw.bad_checks(v)
av = L.all_atom_values(v)
f = L.failing_eqs(av)
print(f"FINAL bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i): v[i] for i in range(L.NVARS)}, open('ordered_state.json', 'w'))
