"""Aim the residue-aware pricer at the candidate clusters ranked by |E| - |inside|.
One real construction per cluster (frame detachment of the bit's two handle variables),
then price.  Anything pricing below 7 is constructed and handed to the checker."""
import sys, os, json, math, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from close4 import *
from pricer import AE, inside, price
W = os.path.dirname(os.path.abspath(__file__)) + '/'
K1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BI = json.load(open(W + 'bitinfo.json'))
TC = json.load(open(W + 'truecost.json'))          # [ |E|, natoms, inside, bit, H1, H2, S ]
def plan(tb, partner):
    d = BI[str(tb)]; q = BI[str(partner)]
    ctrl = {22162: K2, 30213: K1, tb: 1, partner: 1}; det = []
    for pin, T in ((d['xpin'], q['xpin']['C']), (d['ypin'], q['ypin']['C'])):
        X, C, m, H = pin['X'], pin['C'], pin['mult'], pin['H']
        Dv = T - C
        if m == 1: t = 0
        elif math.gcd(P % m, m) == 1: t = (-Dv * pow(P % m, -1, m)) % m
        else: return None
        Xv = T + P * t
        if (Xv - C) % m: return None
        ctrl[X] = Xv; ctrl[H] = (Xv - C) // m; det.append(H)
    return ctrl, det
cands = sorted(TC, key=lambda r: (r[0] - r[2], r[0]))[:int(sys.argv[1]) if len(sys.argv) > 1 else 12]
s1 = [int(b) for b in BI if BI[b]['side'] == 's1']
s2 = [int(b) for b in BI if BI[b]['side'] == 's2']
print('candidate clusters (ranked by |E| - |inside|):', [(r[0], r[2], r[3]) for r in cands], flush=True)
best = (99, None)
for r in cands:
    nE, _, nIn, tb, H1, H2 = r[0], r[1], r[2], r[3], r[4], r[5]
    partner = (s2 if BI[str(tb)]['side'] == 's1' else s1)[0]
    pl = plan(tb, partner)
    if pl is None:
        print('  x_%-6d  plan not constructible (divisibility)' % tb, flush=True); continue
    ctrl, det = pl
    t0 = time.time()
    try: sc, v, nz = closure4(ctrl, detach=det, rounds=12, depth=4)
    except Exception as e:
        print('  x_%-6d  construction failed: %s' % (tb, str(e)[:50]), flush=True); continue
    av = L.all_atom_values(v)
    E = frozenset().union(*[AE[a] for a in (H1, H2)])
    S = inside(E); v0 = [av[a] for a in S]
    pr = price(E, S, v0)
    print('  x_%-6d |E|=%-3d |S|=%-3d resblocked=%-3d maxsat=%-3d PRICE=%-3d  (state %d, %.0fs)'
          % (tb, pr['nE'], pr['nS'], pr['residue_blocked'], pr['max_satisfiable'], pr['price'],
             sc, time.time() - t0), flush=True)
    if pr['price'] < best[0]: best = (pr['price'], tb)
    if pr['price'] < 7:
        out = W + 'PRICED_%d_bit%d.json' % (pr['price'], tb)
        json.dump({f'x_{i}': v[i] for i in range(L.NVARS) if v[i] != 0}, open(out, 'w'))
        print('      *** PRICES BELOW 7 -- written to %s' % out, flush=True)
print('BEST PRICE %d (bit x_%s)' % best, flush=True)
