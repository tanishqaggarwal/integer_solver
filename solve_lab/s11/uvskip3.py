"""From the good closehit2 state: break a41332 [1 eq] and a36244 [4 eqs] to free x_24453 and
   x_3432, hand the MIRROR to them, which releases x_31339 / x_33708 for a33796 / a34580 and
   lets x_33129 / x_37088 close a14445 / a27139 -> a PERFECT matching, cost 5 equations."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
SKIP = {41332, 36244}
sys.set_int_max_str_digits(300000)


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


def drive(v, ctrl, get, tries=6):
    for _ in range(tries):
        r = get(v) % P
        if r == 0:
            return True
        old = v[ctrl]
        v[ctrl] = old + 1
        fwdskip(v)
        s = (get(v) % P - r) % P
        v[ctrl] = old
        fwdskip(v)
        if s == 0:
            return False
        v[ctrl] = (old + (-r) * pow(s, -1, P)) % P
        fwdskip(v)
    return get(v) % P == 0


RES = [('mirror25118', 3432,  lambda v: v[25118]),
       ('mirror3719',  24453, lambda v: v[3719]),
       ('a34580', 33708, lambda v: v[33708] - v[10170]),
       ('a33796', 31339, lambda v: v[31339] - v[6858]),
       ('a14445', 33129, lambda v: v[33129] - v[3757]),
       ('a27139', 37088, lambda v: v[37088] - v[13585])]

v = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'closehit2.json'))).items():
    v[int(k)] = int(x)
fwdskip(v)
print("after breaking the two gates:")
print("   residuals:", [(n, get(v) % P == 0) for n, c, get in RES])
f0 = L.failing_eqs(L.all_atom_values(v))
print(f"   failing={len(f0)} score={L.NEQ-len(f0)}")

# check the two freed controls really move the mirror, and in which pairing
for nm, c, get in RES[:2]:
    old = v[c]
    b = get(v) % P
    v[c] = old + 1
    fwdskip(v)
    print(f"   x{c} moves {nm}: {(get(v)%P - b) % P != 0}")
    v[c] = old
    fwdskip(v)

for rnd in range(8):
    for nm, c, get in RES:
        drive(v, c, get)
    st = [(n, get(v) % P == 0) for n, c, get in RES]
    f = L.failing_eqs(L.all_atom_values(v))
    print(f"  rnd{rnd}: {st} failing={len(f)} score={L.NEQ-len(f)}", flush=True)
    if all(t for _, t in st):
        break

LOCK = {490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742, 22162, 30213,
        8386, 21868, 16441, 28955, 2751, 18751} | {c for _, c, _ in RES}


def bad_now(v):
    return [a for a in range(L.NA) if L.atom_out.get(a) is None and fw.evalpoly(L.polys[a], v) != 0]


best = (len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
for rnd in range(12):
    bad = bad_now(v)
    if not bad:
        break
    prog = False
    for a in sorted(bad, key=lambda a: (len(L.atom2eq.get(a, {})), a)):
        if fw.evalpoly(L.polys[a], v) == 0:
            continue
        cs = [(u, None) for u in L.avars[a]
              if L.definer.get(u) is None and u not in LOCK
              and not any(mm.count(u) > 1 for mm in L.polys[a])]
        cs.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
        try:
            hs, base = deep.handles(v, a, locked=LOCK)
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
    f = L.failing_eqs(L.all_atom_values(v))
    nb = bad_now(v)
    print(f"  close{rnd}: bad={len(nb)} failing={len(f)} score={L.NEQ-len(f)} {nb[:10]}", flush=True)
    if len(f) < best[0]:
        best = (len(f), [x for x in v])
    if not prog or not nb:
        break
print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
json.dump({('x_%d' % i): best[1][i] for i in range(L.NVARS)},
          open(os.path.join(HERE, 'data', 'uvskip4_named.json'), 'w'))
