"""S10 step 48: the last 13 equations on the wire=1 branch.

State: only a37694 (= x_26064 - p, 12 equations) and a39417 (1 equation) are
nonzero -> 39,020.  Find every knob that can move these 13 equations and measure
its collateral, then optimise in equation space.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

ROOTATOM = 37694
BLOCK = {ROOTATOM}
FAILQ = [8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300, 36106,
         36767, 37257, 37666]


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            d = L.definer[u]
            if d in BLOCK:
                continue
            nv = T.solve_lin(d, u, v)
            if nv is not None:
                v[u] = nv
    return v


v = L.load(os.path.join(HERE, 'trade_out.json'))
av = L.all_atom_values(v)
print('nonzero atoms:', [a for a in range(L.NA) if av[a]])
print(f'a37694 = {av[37694]}')
print(f'a39417 = {av[39417]}')
print(f'a39417 src: {L.atom_src[39417][:400]}')
print(f'a37694 src: {L.atom_src[37694]}')

print('\n=== the 13 failing equations, atom by atom ===')
cand_vars = set()
for i in FAILQ:
    m, sq, co = L.eq_atoms[i]
    print(f'\neq {i}  mult={str(m)[:18]} square={sq} atoms={len(co)}')
    for a, c in sorted(co.items()):
        mark = ' <== NONZERO' if av[a] else ''
        out = L.atom_out.get(a)
        print(f'    c={c:<6} a{a:<6} {"GATE" if out else "CHECK":<6} '
              f'neq={len(L.atom2eq.get(a,{})):<3}{mark}  {L.atom_src[a][:80]}')
        for u in L.avars[a]:
            if u not in L.definer:
                cand_vars.add(u)
            else:
                for h in L.avars[L.definer[u]]:
                    if h not in L.definer:
                        cand_vars.add(h)
print(f'\ncandidate free inputs touching these equations: {len(cand_vars)}')

print('\n=== scanning candidates: effect on the 13, and collateral ===', flush=True)
base_fail = set(L.failing_eqs(av))
t0 = time.time()
useful = []
for u in sorted(cand_vars):
    w = list(v); w[u] = w[u] + 1
    fwd_block(w)
    aw = L.all_atom_values(w)
    f2 = set(L.failing_eqs(aw))
    fixed = base_fail - f2
    broke = f2 - base_fail
    if fixed:
        useful.append((len(broke), len(fixed), u, sorted(fixed), sorted(broke)[:8]))
        print(f'  x_{u:<7} fixes {sorted(fixed)} breaks {len(broke)}', flush=True)
useful.sort()
print(f'\nknobs that fix at least one of the 13: {len(useful)}  ({time.time()-t0:.0f}s)')
for b, f, u, fx, bk in useful[:15]:
    print(f'  x_{u:<7} fixes {f} breaks {b}: {fx} / {bk}')
