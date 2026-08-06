"""S10 step 83: do the message bits collapse into a few coefficient GROUPS?

The AD gradient of a29539 showed many distinct boolean inputs carrying the SAME
coefficient (x_91, x_438, x_490, x_33287, x_1203, x_34175 all 3092211522733110...).
If the bits fall into k groups with multiplicities m_1..m_k, then the reachable
residue set is  { sum c_i * delta_i : 0 <= c_i <= m_i },  of size prod(m_i + 1) --
not 2^256.  If that product is small the subset-sum is ENUMERABLE and the target
residue can be tested directly.

For a bit b entering as b*(X - HUGE), the dependence is exactly linear in b, so
the AD coefficient is the exact effect of flipping it.
"""
import os, sys, collections, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
BOOL = set()
for _a, _poly in enumerate(L.polys):
    _ks = list(_poly.items())
    if len(_ks) == 2:
        _sq = [m for m, c in _ks if len(m) == 2 and m[0] == m[1]]
        _li = [m for m, c in _ks if len(m) == 1]
        if _sq and _li and _sq[0][0] == _li[0][0]:
            BOOL.add(_li[0][0])

v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
BFREE = set(u for u in ad.FREE if u in BOOL)
print(f'boolean free inputs: {len(BFREE)}')

CHECKS = [a for a in range(L.NA) if av[a] and a not in atom_out]
print(f'failing checks: {CHECKS}\n')

for c in CHECKS:
    g = ad.grad(c, vm)
    bits = {u: d % P for u, d in g.items() if u in BFREE and d % P}
    if not bits:
        print(f'a{c}: NO boolean input moves it (residue pinned against the message)')
        continue
    groups = collections.Counter(bits.values())
    prod = 1
    for m in groups.values():
        prod *= (m + 1)
        if prod > 10**18:
            break
    print(f'a{c}: {len(bits)} bits move it, {len(groups)} distinct coefficients')
    print(f'   group sizes: {sorted(groups.values(), reverse=True)[:12]}')
    print(f'   reachable-set size prod(m_i+1) = '
          f'{prod if prod <= 10**18 else ">1e18"}')
    tgt = (-av[c]) % P
    print(f'   target shift needed: {str(tgt)[:32]}...')
    if prod <= 5 * 10**6:
        # enumerate the reachable residues exactly
        deltas = sorted(groups.items())
        reach = {0}
        for d, m in deltas:
            nxt = set()
            for r in reach:
                for k in range(m + 1):
                    nxt.add((r + k * d) % P)
            reach = nxt
            if len(reach) > 5 * 10**6:
                break
        print(f'   enumerated {len(reach)} reachable residues; '
              f'TARGET REACHABLE: {tgt in reach}')
        if tgt in reach:
            print('   *** THE MESSAGE CAN HIT THIS CHECK -- reconstruct the bit set!')
    print()

# joint: can the bits hit ALL of them at once?
print('=== joint reachability over the failing checks ===')
mats = {}
for c in CHECKS:
    g = ad.grad(c, vm)
    mats[c] = {u: d % P for u, d in g.items() if u in BFREE and d % P}
allbits = sorted(set().union(*[set(m) for m in mats.values()]))
print(f'bits moving at least one failing check: {len(allbits)}')
sig = collections.Counter(
    tuple(mats[c].get(u, 0) for c in CHECKS) for u in allbits)
print(f'distinct signature vectors among those bits: {len(sig)}')
top = sig.most_common(8)
for s, n in top:
    print(f'   x{n} bits share signature {[str(x)[:12] for x in s]}')
prod = 1
for _, n in sig.items():
    prod *= (n + 1)
    if prod > 10**18:
        prod = None; break
print(f'joint reachable-set size prod(m_i+1) = {prod if prod else ">1e18"}')
