"""Recover the selector tree STRUCTURALLY from my own frame, with no dependence on other
agents' chains.

Every defined variable v has a support bitmask over Frame's free inputs (`fr.sup`).  Restrict
that support to the 256 selector variables (the keys of leafpins.json, all of which are pure
free inputs).  The distinct nonempty restricted supports form a LAMINAR family — that family is
the OR-tree, recovered without parsing a single gate.

Output: runs/seltree.json = {roots, blocks (sorted by size), parent, depth, leaf order}.
"""
import os, sys, json, time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import optN
from optN import fr, FR0

SEL = sorted(int(k) for k in json.load(open(os.path.join(HERE, 'leafpins.json'))))
SELSET = set(SEL)
sidx = {s: i for i, s in enumerate(SEL)}
NS = len(SEL)

# selector-only mask for each free input
fidx = fr.fidx
selmask_free = {}
for s in SEL:
    assert s in fidx, ('selector not free in frame', s)
    selmask_free[fidx[s]] = 1 << sidx[s]


def sel_support(mask):
    """map a Frame free-input bitmask -> bitmask over the 256 selectors"""
    out = 0
    m = mask
    while m:
        b = m & -m
        j = b.bit_length() - 1
        q = selmask_free.get(j)
        if q:
            out |= q
        m ^= b
    return out


def popcount(x):
    return bin(x).count('1')


def main():
    t0 = time.time()
    blocks = defaultdict(list)   # selector-support bitmask -> list of vars carrying it
    for v in fr.order:
        m = sel_support(fr.sup[v])
        if m:
            blocks[m].append(v)
    for a in fr.checks:
        m = sel_support(fr.csup[a])
        if m:
            blocks[m].append(('a', a))
    print('distinct nonempty selector-supports: %d   (%.1fs)' % (len(blocks), time.time() - t0))

    keys = sorted(blocks, key=lambda m: (popcount(m), m))
    sizes = defaultdict(int)
    for m in keys:
        sizes[popcount(m)] += 1
    print('by size:', dict(sorted(sizes.items())[:20]), '...' if len(sizes) > 20 else '')
    print('largest sizes:', sorted(sizes.items())[-8:])

    # laminarity test over the distinct supports
    K = keys
    bad = 0
    bysize = sorted(K, key=popcount)
    for i, a in enumerate(bysize):
        for b in bysize[i + 1:]:
            if popcount(b) > popcount(a):
                if (a & b) and (a & b) != a:
                    bad += 1
                    if bad <= 3:
                        print('NON-LAMINAR pair popcounts', popcount(a), popcount(b),
                              'inter', popcount(a & b))
    print('non-laminar pairs:', bad, 'of', len(K) * (len(K) - 1) // 2)

    # build the tree over the laminar family
    parent = {}
    for i, a in enumerate(bysize):
        best = None
        for b in bysize:
            if b == a or popcount(b) <= popcount(a):
                continue
            if (a & b) == a:
                if best is None or popcount(b) < popcount(best):
                    best = b
        parent[a] = best
    roots = [a for a in bysize if parent[a] is None]
    print('ROOTS:', [(popcount(r), len(blocks[r])) for r in roots])

    depth = {}

    def dep(m):
        if m in depth:
            return depth[m]
        p = parent[m]
        depth[m] = 0 if p is None else dep(p) + 1
        return depth[m]

    for m in bysize:
        dep(m)
    bd = defaultdict(int)
    for m in bysize:
        bd[depth[m]] += 1
    print('blocks by depth:', dict(sorted(bd.items())))

    singletons = [m for m in bysize if popcount(m) == 1]
    print('singleton blocks (individual selectors seen alone):', len(singletons))

    out = dict(
        selectors=SEL,
        roots=[[popcount(r), [SEL[i] for i in range(NS) if r >> i & 1]] for r in roots],
        blocks=[dict(size=popcount(m), depth=depth[m],
                     members=[SEL[i] for i in range(NS) if m >> i & 1],
                     nvars=len(blocks[m]))
                for m in bysize],
        nonlaminar=bad,
    )
    json.dump(out, open(os.path.join(HERE, 'runs', 'seltree.json'), 'w'))
    print('wrote runs/seltree.json  (%.1fs)' % (time.time() - t0))


if __name__ == '__main__':
    main()
