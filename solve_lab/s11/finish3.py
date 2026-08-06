"""Close the last three: a14445, a27139 directly; a26719 by a p-multiple shift making
   8640431 | (x_12000 / p)."""
import sys, os, json, time, random, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, uv01, uv01full
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
M = 8640431
print("8640431 factors:", [(d, M // d) for d in range(2, 4000) if M % d == 0][:6])

v = [0] * L.NVARS
d = json.load(open(os.path.join(HERE, 'data', 'three.json')))
for k, x in d.items():
    v[int(k)] = int(x)
fw.forward(v)
print("start bad:", fw.bad_checks(v))

# 1) direct closes
for a, t in [(14445, 33129), (27139, 37088)]:
    x = fw.solve_lin(a, t, v)
    if x is not None:
        v[t] = x
        fw.forward(v)
    print(f"  closed a{a} via x{t}: {fw.evalpoly(L.polys[a],v)==0}")
print("bad now:", fw.bad_checks(v))
print(f"  x3719%P==0:{v[3719]%P==0} x25118%P==0:{v[25118]%P==0}")

# 2) the 8640431 condition
alpha = v[3719] // P
beta = v[25118] // P
gam = v[12000] // P
print(f"  alpha%M={alpha%M} beta%M={beta%M} gamma%M={gam%M}  (need gamma%M==0)")
print(f"  check x12000 == 9974121*x3719 + 15683097*x25118 :",
      v[12000] == 9974121 * v[3719] + 15683097 * v[25118])

# measure response of (x3719, x25118) to a +P shift of each mirror control
CTRL = uv01full.MIRROR_CTRL
resp = {}
for c in CTRL:
    old = v[c]
    v[c] = old + P
    fw.forward(v)
    ok = (v[3719] % P == 0 and v[25118] % P == 0)
    da = v[3719] // P - alpha
    db = v[25118] // P - beta
    dg = v[12000] // P - gam
    resp[c] = (da % M, db % M, dg % M, ok)
    v[c] = old
    fw.forward(v)
    print(f"  x{c} += P : d(alpha)={da%M} d(beta)={db%M} d(gamma)={dg%M} modp-mirror-kept={ok}")

# solve  gamma + sum_c k_c * dg_c  ==  0  (mod M)
gs = [resp[c][2] for c in CTRL]
g = 0
for x in gs:
    g = math.gcd(g, x)
g = math.gcd(g, M)
print(f"  gcd of dgamma with M = {g};  (-gamma) % g == 0 ? {(-gam) % M % g == 0 if g else None}")
sol = None
for i, c in enumerate(CTRL):
    dgc = gs[i]
    gg = math.gcd(dgc, M)
    if gg and (-gam) % M % gg == 0:
        mm = M // gg
        k = ((-gam) % M // gg) * pow((dgc // gg) % mm, -1, mm) % mm
        sol = (c, k)
        break
print("  single-control solution:", sol)
if sol:
    c, k = sol
    v[c] = v[c] + k * P
    fw.forward(v)
    print(f"  applied x{c} += {k}*P")
    print(f"    x3719%P==0:{v[3719]%P==0} x25118%P==0:{v[25118]%P==0} gamma%M={(v[12000]//P)%M}")
    bad = fw.bad_checks(v)
    print("  bad:", bad)
    for a in bad:
        uv01full.close_one(v, a, set())
    bad = fw.bad_checks(v)
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"FINAL bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad}")
    json.dump({str(i): v[i] for i in range(L.NVARS)}, open(os.path.join(HERE, 'data', 'finish3.json'), 'w'))
