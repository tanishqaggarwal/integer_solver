"""Full constructive assembly: structural solve + deep cone-handle closure."""
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
TH = {int(k): x for k, x in json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'theta_solveB.json')).items()}

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
b = fw.bad_checks(v)
print("bad after structural:", len(b), b)
print("  mod-p nonzero:", [a for a in b if fw.evalpoly(L.polys[a], v) % P != 0])


def shallow(a, locked):
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked:
            continue
        if any(mm.count(u) > 1 for mm in L.polys[a]):
            continue
        out.append(u)
    out.sort(key=lambda u: (NAT[u], u))
    return out


t0 = time.time()
for rnd in range(40):
    bad = fw.bad_checks(v)
    if not bad:
        break
    prog = False
    for a in bad:
        if fw.evalpoly(L.polys[a], v) == 0:
            continue
        done = False
        # 1) shallow: direct free var in the atom
        for t in shallow(a, LOCK):
            x = fw.solve_lin(a, t, v)
            if x is not None and x != v[t]:
                old = v[t]
                v[t] = x
                fw.forward(v)
                if fw.evalpoly(L.polys[a], v) == 0:
                    done = prog = True
                    break
                v[t] = old
                fw.forward(v)
        if done:
            continue
        # 2) deep: exact linear handle anywhere in the atom's cone
        try:
            hs, base = deep.handles(v, a, locked=LOCK)
        except Exception:
            continue
        for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0])):
            if d == 0 or base % d:
                continue
            old = v[t]
            v[t] = old - base // d
            fw.forward(v)
            if fw.evalpoly(L.polys[a], v) == 0:
                prog = True
                done = True
                break
            v[t] = old
            fw.forward(v)
    bad2 = fw.bad_checks(v)
    print(f"  round{rnd}: bad={len(bad2)} ({time.time()-t0:.0f}s) {bad2[:14]}", flush=True)
    if not prog:
        break
fw.forward(v)
b = fw.bad_checks(v)
av = L.all_atom_values(v)
f = L.failing_eqs(av)
print(f"FINAL bad_checks={len(b)} failing_eqs={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i): v[i] for i in range(L.NVARS)}, open('assemble.json', 'w'))
