"""Two-parameter p-shift:  x31339 += k*P , x33708 += l*P.
   Both keep x3719 == x25118 == 0 (mod P).  gamma(k,l) = x_12000/P is a quadratic FORM.
   Solve gamma == 0 mod 8640431 = 53 * 163027 by CRT."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from quad8640431 import quad_roots, load
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
M = 8640431
QS = [53, 163027]


def gamma_form(v, c1=31339, c2=33708):
    """interpolate gamma(k,l) = a00 + a10 k + a01 l + a20 k^2 + a11 k l + a02 l^2 (mod M)"""
    b1, b2 = v[c1], v[c2]
    pts = [(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1)]
    ys = []
    for k, l in pts:
        v[c1] = b1 + k * P
        v[c2] = b2 + l * P
        fw.forward(v)
        if v[3719] % P or v[25118] % P:
            v[c1], v[c2] = b1, b2
            fw.forward(v)
            return None
        ys.append((v[12000] // P) % M)
    v[c1], v[c2] = b1, b2
    fw.forward(v)
    a00 = ys[0]
    # k-direction
    a20 = ((ys[2] - 2 * ys[1] + ys[0]) * pow(2, -1, M)) % M
    a10 = (ys[1] - ys[0] - a20) % M
    a02 = ((ys[4] - 2 * ys[3] + ys[0]) * pow(2, -1, M)) % M
    a01 = (ys[3] - ys[0] - a02) % M
    a11 = (ys[5] - a00 - a10 - a01 - a20 - a02) % M
    return (a00, a10, a01, a20, a11, a02)


def evalf(F, k, l, m):
    a00, a10, a01, a20, a11, a02 = F
    return (a00 + a10 * k + a01 * l + a20 * k * k + a11 * k * l + a02 * l * l) % m


def solve_mod(F, q, limit=None):
    """all (k,l) mod q with gamma==0; iterate k, solve quadratic in l"""
    a00, a10, a01, a20, a11, a02 = [x % q for x in F]
    out = []
    rng = range(q)
    for k in rng:
        A = a02
        B = (a01 + a11 * k) % q
        C = (a00 + a10 * k + a20 * k * k) % q
        for l in quad_roots(A, B, C, q):
            out.append((k, l))
            if limit and len(out) >= limit:
                return out
    return out


if __name__ == '__main__':
    v = load('three.json')
    print("state bad:", fw.bad_checks(v), "failing:", len(L.failing_eqs(L.all_atom_values(v))))
    F = gamma_form(v)
    print("gamma form:", F)
    if F is None:
        print("shift breaks the mirror"); sys.exit()
    s53 = solve_mod(F, 53)
    print(f"solutions mod 53: {len(s53)}")
    s163 = solve_mod(F, 163027, limit=40)
    print(f"solutions mod 163027 (first 40): {len(s163)}")
    if not s53 or not s163:
        print("no CRT solution from this mirror branch"); sys.exit()
    inv53 = pow(163027, -1, 53)
    inv163 = pow(53, -1, 163027)
    b1, b2 = v[31339], v[33708]
    tried = 0
    for (k1, l1) in s53[:8]:
        for (k2, l2) in s163[:8]:
            k = (k1 * 163027 * inv53 + k2 * 53 * inv163) % M
            l = (l1 * 163027 * inv53 + l2 * 53 * inv163) % M
            assert evalf(F, k, l, M) == 0
            v[31339] = b1 + k * P
            v[33708] = b2 + l * P
            fw.forward(v)
            tried += 1
            gam = (v[12000] // P) % M
            bad = fw.bad_checks(v)
            f = L.failing_eqs(L.all_atom_values(v))
            print(f"  k,l -> gamma%M={gam} mirror0={v[3719]%P==0 and v[25118]%P==0} "
                  f"bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:10]}", flush=True)
            if gam == 0 and len(bad) < 3:
                json.dump({str(i): v[i] for i in range(L.NVARS)},
                          open(os.path.join(HERE, 'data', 'quad2_hit.json'), 'w'))
                print("   SAVED")
            if tried >= 12:
                break
        if tried >= 12:
            break
    v[31339], v[33708] = b1, b2
    fw.forward(v)
