import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
BITS = (542, 47, 438, 91)
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
TH = {int(k): v for k, v in json.load(open('theta_solveA.json')).items()}

v = [0] * L.NVARS
for b in BITS:
    v[b] = 1
for k, x in TH.items():
    v[k] = x
fw.forward(v)
print("targets:", [v[3719] % P, v[25118] % P, v[25614] % P, v[34220] % P,
                   (v[12186]-v[1308]) % P, (v[24908]-v[19083]) % P])
print("x15298=%d x5647=%d x34606=%d" % (v[15298], v[5647], v[34606]))

# kill first core exactly
v[14853] = v[12186]
v[16742] = v[24908]
fw.forward(v)
print("n=%d m=%d" % (v[29322], v[3558]))

# arithmetic cluster: x37892 = x30213, x13682 = x22162 in channel U=V=1
v[30213] = G0
v[22820] = 0
v[7497] = (C0 + MM * G0) // P
v[22162] = -C0B
v[14393] = 0
v[11436] = 0
fw.forward(v)
print("x37892==x30213:", v[37892] == v[30213], " x13682==x22162:", v[13682] == v[22162])
print("a688=%s a1618=%s a40608=%s" % (fw.evalpoly(L.polys[688], v), fw.evalpoly(L.polys[1618], v),
                                      fw.evalpoly(L.polys[40608], v)))
b = fw.bad_checks(v)
print("bad before repair:", len(b), b[:30])
print("  mod-p nonzero:", [a for a in b if fw.evalpoly(L.polys[a], v) % P != 0])

LOCK = set(BITS) | set(TH) | {14853, 16742, 30213, 22820, 7497, 22162, 14393, 11436}


def free_cands(a, locked):
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
for it in range(60):
    b = fw.bad_checks(v)
    if not b:
        break
    prog = False
    for a in b:
        if fw.evalpoly(L.polys[a], v) == 0:
            continue
        for t in free_cands(a, LOCK):
            x = fw.solve_lin(a, t, v)
            if x is not None and x != v[t]:
                old = v[t]
                v[t] = x
                fw.forward(v)
                if fw.evalpoly(L.polys[a], v) == 0:
                    prog = True
                    break
                v[t] = old
                fw.forward(v)
    if not prog:
        break
fw.forward(v)
b = fw.bad_checks(v)
av = L.all_atom_values(v)
f = L.failing_eqs(av)
print(f"AFTER REPAIR: bad_checks={len(b)} failing_eqs={len(f)} score={L.NEQ-len(f)} ({time.time()-t0:.0f}s)")
print("bad:", b)
print("  mod-p nonzero:", [a for a in b if fw.evalpoly(L.polys[a], v) % P != 0])
json.dump({str(i): v[i] for i in range(L.NVARS)}, open('final_state.json', 'w'))
