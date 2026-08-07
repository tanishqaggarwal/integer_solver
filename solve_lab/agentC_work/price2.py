"""Corrected: price the cluster that ACTUALLY carries the residual in each construction.
S = the atoms that are genuinely nonzero in the state (the defect support);
E = the union of their equations; v0 = their measured values.
Sanity gate: price must equal the observed failing count, else the cluster is not closed
(the defect leaks into equations outside E) and the number is not trustworthy."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from close4 import *
from pricer import AE, inside, price
W = os.path.dirname(os.path.abspath(__file__)) + '/'
K1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BI = json.load(open(W + 'bitinfo.json'))
TC = json.load(open(W + 'truecost.json'))
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
# --- re-calibrate the corrected form on the deliverable
d = json.load(open(W + '../best/new_instance_partial_39026.json'))
v = [0] * L.NVARS
for k, val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
av = L.all_atom_values(v); fails = L.failing_eqs(av)
S = [a for a in range(L.NA) if av[a]]
E = frozenset().union(*[AE[a] for a in S])
Sc = inside(E)
pr = price(E, Sc, [av[a] for a in Sc])
print('CALIBRATION (corrected form): |E|=%d |S|=%d price=%d   observed failing=%d  -> %s'
      % (pr['nE'], pr['nS'], pr['price'], len(fails), 'MATCH' if pr['price'] == len(fails) else 'MISMATCH'), flush=True)
s1 = [int(b) for b in BI if BI[b]['side'] == 's1']
s2 = [int(b) for b in BI if BI[b]['side'] == 's2']
cands = sorted(TC, key=lambda r: (r[0] - r[2], r[0]))[:int(sys.argv[1]) if len(sys.argv) > 1 else 8]
print('\nbit    |E|  |S| resblk maxsat PRICE  observed  closed?', flush=True)
for r in cands:
    tb = r[3]
    partner = (s2 if BI[str(tb)]['side'] == 's1' else s1)[0]
    pl = plan(tb, partner)
    if pl is None: continue
    ctrl, det = pl
    try: sc, vv, nz = closure4(ctrl, detach=det, rounds=12, depth=4)
    except Exception: continue
    av2 = L.all_atom_values(vv); f2 = L.failing_eqs(av2)
    S2 = [a for a in range(L.NA) if av2[a]]
    if not S2: continue
    E2 = frozenset().union(*[AE[a] for a in S2]); Sc2 = inside(E2)
    p2 = price(E2, Sc2, [av2[a] for a in Sc2])
    closed = (p2['price'] == len(f2))
    print('x_%-6d %-4d %-4d %-6d %-6d %-6d %-9d %s'
          % (tb, p2['nE'], p2['nS'], p2['residue_blocked'], p2['max_satisfiable'],
             p2['price'], len(f2), 'yes' if closed else 'NO (leaks)'), flush=True)
    if closed and p2['price'] < 7:
        out = W + 'PRICED_%d_bit%d.json' % (p2['price'], tb)
        json.dump({f'x_{i}': vv[i] for i in range(L.NVARS) if vv[i] != 0}, open(out, 'w'))
        print('   *** below 7 and closed -> %s' % out, flush=True)
