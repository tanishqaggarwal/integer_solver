"""S11 step 98: the exact map from EVERY free input to the six values.

Part XXVII said all seven quantities in A and B are literal constants.  That was
too strong, and the support computation says why:

    w1 = x12186   computed, 417-input support        w2 = x16742   FREE
    w3 = x14853   FREE                               w4 = x24908   computed, 99 inputs
    w5 = x22162   FREE                               w6 = x30213   FREE
    K  = x24453   constant (empty support)           <- only this one is really fixed

and perturbing x22152, x33462, x6418 or x12553 moves NONE of w1, w2, w3, w4.  Those
four literals are separate variables that merely happen to carry the same residues at
the current state; they do not feed A or B.  So the values are constrained
by congruences to COMPUTED targets, not pinned to constants, and the question becomes
what actually moves them.

One forward-AD pass per free input gives the exact derivative of every value --
7,273 passes, affordable -- so the linear map from the whole free-input space to the
six-dimensional w space can be written down exactly.

Usage: valjac.py START END [state.json]
"""
import os, sys, time, json, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from chunk import sweep, load
P = ad.P
src = sys.argv[3] if len(sys.argv) > 3 else 'PIN_39013.json'
tag = 'coord_' + os.path.basename(src).replace('.json', '')
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
vm = [x % P for x in v]
KEY = [12186, 16742, 14853, 24908, 22162, 30213, 19083, 1308, 35389, 6671]
NAMES = dict(zip(KEY, ['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'y1_tgt', 'x2_tgt',
                       'A', 'B']))
FREE = [t for t in range(L.NVARS) if t not in L.definer]
print('%s: %d free inputs; tracking %s'
      % (src, len(FREE), [NAMES[k] for k in KEY]), flush=True)


def dcol(u):
    """d(variable)/d(x_u) mod p for the tracked variables, one forward pass."""
    dv = {u: 1}
    for t in ad.ORDER:
        a = L.definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0:
            dv[t] = 0
            continue
        s = 0
        for w in L.avars[a]:
            if w == t:
                continue
            dw = dv.get(w, 0)
            if dw:
                s += ad.dpart(a, w, vm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    return {k: dv.get(k, 0) % P for k in KEY}


def evaluate(u):
    d = dcol(u)
    return {'u': u, 'd': {str(k): str(x) for k, x in d.items() if x},
            'n': sum(1 for x in d.values() if x)}


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(FREE)
sweep(tag, FREE, evaluate, start, min(end, len(FREE)), keyfn=str, budget=540)
rs = load(tag)
print('\n%d free inputs measured' % len(rs), flush=True)
hits = collections.Counter()
movers = collections.defaultdict(list)
for r in rs:
    for k, x in r['d'].items():
        hits[NAMES[int(k)]] += 1
        movers[NAMES[int(k)]].append(r['u'])
print('free inputs that move each value:')
for nm in ['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'y1_tgt', 'x2_tgt', 'A', 'B']:
    print('   %-7s %-5d   sample %s' % (nm, hits[nm], movers[nm][:10]))
json.dump({k: v2 for k, v2 in movers.items()},
          open(os.path.join(HERE, 'valmovers.json'), 'w'))
print('saved valmovers.json')
