"""S11 step 96: coherent ONE-HOT re-selection, not single-bit flips.

Reading the atoms of a boolean free input shows what the bits actually are:

    a20545  2*b - 2*b*b                        the boolean constraint
    a24804  x18232 = b        a24805  x36695 = 1 - b      the bit and its complement
    a21775  b*x5803 - CONST*b - 12107359*x22874           IF b THEN x5803 = CONST
    a35126  b*x38738 - CONST2*b - x12204                  IF b THEN x38738 = CONST2

so each bit is a **selector over a table of constants** -- a multiplexer, exactly what
a windowed scalar multiplication uses to pick a precomputed point.  That explains the
census: flipping ONE bit breaks the one-hot invariant and costs 30-56, and only two
distinct (dA, dB) outcomes exist among 400 costing bits, because turning any bit on
fires the same OR.

The coherent move is a SWAP: turn the currently-selected bit off and another bit of
the same group on.  That changes which table entry is read -- which changes the point
being added -- while keeping the one-hot invariant.  No single-bit census can see it,
and nothing in this lab has tried it.

Groups are recovered from the gated pins: two bits are in the same group when their
pins select the same wire.

Usage: onehot.py START END [state.json]
"""
import os, sys, collections, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from chunk import sweep, load
P = ad.P
src = sys.argv[3] if len(sys.argv) > 3 else 'PIN_39013.json'
tag = 'onehot_' + os.path.basename(src).replace('.json', '')
base = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(base, rounds=6)
BASE = L.NEQ - len(L.failing_eqs(L.all_atom_values(base)))
A0, B0 = base[35389] % P, base[6671] % P
FREE = [t for t in range(L.NVARS) if t not in L.definer]
BITS = [t for t in FREE if base[t] in (0, 1)]
print('%s: score %d; %d boolean free inputs' % (src, BASE, len(BITS)), flush=True)

# a gated pin on bit b:  b*x - C*b - handle    ->  b selects wire x
sel = collections.defaultdict(set)          # wire -> bits that select it
pins = collections.defaultdict(dict)        # bit -> {wire: constant}
for b in BITS:
    for a in L.var_atoms[b]:
        if a in L.atom_out:
            continue
        pl = L.polys[a]
        if len(pl) > 4:
            continue
        for m, c in pl.items():
            if len(m) == 2 and b in m:
                w = m[0] if m[1] == b else m[1]
                if w == b:
                    continue
                lin = [(mm, cc) for mm, cc in pl.items()
                       if len(mm) == 1 and mm[0] == b]
                if lin:
                    sel[w].add(b)
                    pins[b][w] = (-lin[0][1]) * pow(c, -1, P) % P
GROUPS = {w: sorted(bs) for w, bs in sel.items() if len(bs) > 1}
print('%d selector groups (wire -> bits); sizes %s'
      % (len(GROUPS), sorted(collections.Counter(len(v) for v in GROUPS.values())
                             .items())[:8]), flush=True)
ON = [b for b in BITS if base[b] == 1]
print('%d bits are currently 1' % len(ON), flush=True)

CAND = []
for w, bs in sorted(GROUPS.items()):
    on = [b for b in bs if base[b] == 1]
    off = [b for b in bs if base[b] == 0]
    if not on:
        continue
    for a in on:
        for b in off:
            CAND.append((a, b))
print('%d coherent swaps available' % len(CAND), flush=True)


def evaluate(spec):
    a, b = spec
    v = list(base)
    v[a] = 0
    v[b] = 1
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    out = {'off': a, 'on': b, 'score': s,
           'nz': sum(1 for x in range(L.NA) if x not in L.atom_out and av[x]),
           'A0': int(v[35389] % P == 0), 'B0': int(v[6671] % P == 0),
           'movedA': int(v[35389] % P != A0), 'movedB': int(v[6671] % P != B0)}
    if s > BASE:
        T.save(v, os.path.join(HERE, 'OH_%d_%d_%d.json' % (s, a, b)))
    return out


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(CAND)
sweep(tag, CAND, evaluate, start, min(end, len(CAND)),
      keyfn=lambda s: '%d_%d' % s, budget=540)
rs = load(tag)
if rs:
    rs.sort(key=lambda r: -r['score'])
    print('\ntop swaps from %d:' % BASE)
    for r in rs[:15]:
        print('   x%-6d off / x%-6d on -> %-6d  checks %-3s  A moved %s  B moved %s'
              % (r['off'], r['on'], r['score'], r.get('nz'), r.get('movedA'),
                 r.get('movedB')))
    mv = [r for r in rs if r.get('movedA') or r.get('movedB')]
    print('\nswaps that move A or B: %d;  that zero one of them: %d'
          % (len(mv), sum(1 for r in rs if r.get('A0') or r.get('B0'))))
