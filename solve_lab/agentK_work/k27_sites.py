#!/usr/bin/env python3
"""K27: equation footprint of every candidate 2-parameter injection site.

A full solve needs fold(S) == target.  Failing that, the cheapest partial is obtained by
breaking the *smallest* set of equations while buying two free scalars, then inverting the
target down to that site (the stage law inverts in closed form).  This measures the cost of
every site: leaf pins, stage input-slot pins, and every other residual atom pair."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P

C = CascadeP()
atomeq = collections.defaultdict(set)
for e, row in enumerate(C.eqrows):
    for k, j in row: atomeq[j].add(e)
cnt = collections.Counter(len(v) for v in atomeq.values())
print('equations-per-atom histogram:', sorted(cnt.items())[:12], '...max', max(cnt))
solo = [j for j, v in atomeq.items() if len(v) == 1]
print('atoms appearing in exactly 1 equation:', len(solo))

D = json.load(open(K + '/points.json'))
leafsel = {l['sel']: (l['wx'], l['wy']) for l in D['leaves']}
byvar = C.var2atoms

# a site is a pair of atoms whose breaking frees exactly the x and y of one wire pair
def cost(js):
    u = set()
    for j in js: u |= atomeq[j]
    return len(u)

print()
print('--- leaf pin sites ---')
best = []
pinpat = re.compile(r'^\(\(x(\d+)\*\(x(\d+)-\d{20,}\)\)-')
leafpin = {}
for j, nm in enumerate(C.names):
    m = pinpat.match(nm)
    if m: leafpin.setdefault(int(m.group(1)), {})[int(m.group(2))] = j
for s, wm in leafpin.items():
    wx, wy = leafsel[s]
    js = [wm[wx], wm[wy]]
    best.append((cost(js), 'leafpin', s, js))
best.sort()
for b in best[:8]: print('   cost %2d  leaf sel x%d  eqs %s' % (b[0], b[2], sorted(set().union(*[atomeq[j] for j in b[3]]))))

print()
print('--- all 2-atom sites over residual atoms on a common wire pair is too broad;')
print('    instead: cheapest individual residual atoms ---')
resset = set(C.aidx[a] for a in C.E.res)
cheap = sorted(((len(atomeq[j]), j) for j in resset))[:25]
for c, j in cheap: print('   %2d eqs  %s' % (c, C.names[j][:100]))
json.dump({'leafpin_costs': [[b[0], b[2]] for b in best]}, open(K + '/sitecost.json', 'w'))
print()
print('min leaf-pin site cost:', best[0][0], ' (baseline deliverable costs 7)')
