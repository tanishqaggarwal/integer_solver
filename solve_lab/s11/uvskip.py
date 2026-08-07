"""Beat the deficit by BREAKING two cheap gate atoms instead of two expensive checks.
     a41332 [1 eq]  -> frees x_24453
     a36244 [4 eqs] -> frees x_3432
   Both move the mirror, which had no non-bit control.  Total cost 5 equations.
"""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from gfp import gauss_solve
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
LD = json.load(open(os.path.join(HERE, 'data', 'loads.json')))['loads']
C0 = L.polys[688][()]
MM = 8863713
G0 = (-C0 * pow(MM, -1, P)) % P
C0B = L.polys[1618][()]
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
SKIP = {41332, 36244}
MCTRL = [24453, 3432]
LINKS = [(21050, 16441, lambda v: v[16441] - v[4920]),
         (40065, 28955, lambda v: v[28955] - v[11408]),
         (7881,  2751,  lambda v: v[2751] - v[1085]),
         (26839, 18751, lambda v: v[18751] - v[33091]),
         (14445, 33129, lambda v: v[33129] - v[3757]),
         (34580, 33708, lambda v: v[33708] - v[10170]),
         (27139, 37088, lambda v: v[37088] - v[13585]),
         (33796, 31339, lambda v: v[31339] - v[6858])]


def fwdskip(v, skip=SKIP):
    for comp in fw.ORDER:
        if len(comp) == 1:
            u = comp[0]
            if L.definer[u] in skip:
                continue
            x = fw.solve_lin(L.definer[u], u, v)
            if x is not None:
                v[u] = x
        else:
            for _ in range(40):
                ch = False
                for u in comp:
                    if L.definer[u] in skip:
                        continue
                    x = fw.solve_lin(L.definer[u], u, v)
                    if x is not None and x != v[u]:
                        v[u] = x
                        ch = True
                if not ch:
                    break
    return v


def drive(v, ctrl, get, target=0, tries=5):
    for _ in range(tries):
        r = (get(v) - target) % P
        if r == 0:
            return True
        old = v[ctrl]
        v[ctrl] = old + 1
        fwdskip(v)
        s = ((get(v) - target) % P - r) % P
        v[ctrl] = old
        fwdskip(v)
        if s == 0:
            return False
        v[ctrl] = (old + (-r) * pow(s, -1, P)) % P
        fwdskip(v)
    return (get(v) - target) % P == 0


def arithmetic(v):
    v[22162] = 0
    v[30213] = 0
    v[8386] = 0
    v[21868] = 0
    fwdskip(v)
    drive(v, 19750, lambda vv: vv[37892], G0)
    num = C0 + MM * v[37892]
    v[7497] = num // P if num % P == 0 else 0
    v[22820] = 0
    fwdskip(v)
    v[14853] = (-C0B) % P
    fwdskip(v)
    dd = v[13682] + C0B
    v[14393] = 0
    v[11436] = (dd // P) if dd % P == 0 else 0
    fwdskip(v)
    drive(v, 14515, lambda vv: vv[1308] - vv[14853], 0)
    v[16742] = v[19083]
    fwdskip(v)


def mirror(v, iters=6):
    """TRIANGULAR: x_25118 is independent of x_24453 (which enters only x_23776),
       and x_3719 is LINEAR in x_24453 with coefficient x_3090^2."""
    for it in range(iters):
        drive(v, 3432, lambda vv: vv[25118], 0)
        drive(v, 24453, lambda vv: vv[3719], 0)
        if not (v[3719] % P or v[25118] % P):
            return True
    return not (v[3719] % P or v[25118] % P)


def bad_checks_skip(v):
    return [a for a in range(L.NA)
            if (L.atom_out.get(a) is None or a in SKIP) and fw.evalpoly(L.polys[a], v) != 0]


def close_all(v, locked, rounds=15, verbose=True):
    for rnd in range(rounds):
        bad = [a for a in bad_checks_skip(v) if a not in SKIP]
        if not bad:
            break
        prog = False
        for a in sorted(bad, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
            if fw.evalpoly(L.polys[a], v) == 0:
                continue
            cs = []
            for u in L.avars[a]:
                if L.definer.get(u) is None and u not in locked and \
                        not any(mm.count(u) > 1 for mm in L.polys[a]):
                    cs.append((u, None))
            cs.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
            try:
                hs, base = deep.handles(v, a, locked=locked)
                cs += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
            except Exception:
                pass
            for t, d in cs:
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
                fwdskip(v)
                if fw.evalpoly(L.polys[a], v) == 0:
                    prog = True
                    break
                v[t] = old
                fwdskip(v)
        nb = [a for a in bad_checks_skip(v) if a not in SKIP]
        f = L.failing_eqs(L.all_atom_values(v))
        if verbose:
            print(f"    close{rnd}: bad={len(nb)} failing={len(f)} score={L.NEQ-len(f)} {nb[:10]}", flush=True)
        if not prog or set(nb) == set(bad):
            break
    return v


if __name__ == '__main__':
    eqs = set()
    for a in SKIP:
        eqs |= set(L.atom2eq.get(a, {}))
    print(f"broken gates {sorted(SKIP)} occupy {len(eqs)} equations: {sorted(eqs)}")
    v = [0] * L.NVARS
    for b in (490, 91):
        v[b] = 1
    fwdskip(v)
    for b in (490, 91):
        for rec in LD.get(str(b), []):
            a, X = rec[0], rec[1]
            x = fw.solve_lin(a, X, v)
            if x is not None:
                v[X] = x
    fwdskip(v)
    arithmetic(v)
    for a, c, g in LINKS:
        drive(v, c, g, 0)
    m = mirror(v)
    arithmetic(v)
    for a, c, g in LINKS:
        drive(v, c, g, 0)
    print(f"mirror solved: {m}  x3719%p==0:{v[3719]%P==0} x25118%p==0:{v[25118]%P==0}")
    LOCK = {490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742, 22162, 30213,
            8386, 21868} | set(MCTRL) | {c for _, c, _ in LINKS}
    close_all(v, LOCK)
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"FINAL failing={len(f)} score={L.NEQ-len(f)}")
    sys.set_int_max_str_digits(300000)
    json.dump({('x_%d' % i): v[i] for i in range(L.NVARS)},
              open(os.path.join(HERE, 'data', 'uvskip_named.json'), 'w'))
