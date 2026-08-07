#!/usr/bin/env python3
"""K34: localize the B-half gap exactly.

For every wire whose leaf support the descent recovered, the model predicts
   value(w) = composition of the live leaves in support(w).
Compare that against what the closure actually put on the wire, and report the SHALLOWEST
wire that disagrees while all of its support's sub-wires agree.  That is the stage where the
model and the equations part company."""
import sys, os, json, re
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import mux as MUX
import fold as FD
from k26_drive import drive, rootpair, C, P

ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow, selof = {}, {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
exp2sel = {e: s for s, e in sel2exp.items()}
S = FD.SHIFT

leafsel, leafcoord = {}, {}
for l in D['leaves']:
    leafsel[l['wx']] = l['sel']; leafsel[l['wy']] = l['sel']
    leafcoord[l['wx']] = 'x'; leafcoord[l['wy']] = 'y'
gatedpat = re.compile(r'^\(x(\d+)\*x(\d+)\)$')
memo = {}


def support(w, depth=0):
    if w in memo: return memo[w]
    if w in leafsel:
        memo[w] = frozenset({leafsel[w]}); return memo[w]
    memo[w] = frozenset()
    if depth > 40: return memo[w]
    out = set()
    for z, coef in MUX.source_of(w):
        if z == 'CONST': continue
        for kind, t in MUX.mux_terms(z):
            if kind == 'gated':
                m = gatedpat.match(t)
                if m:
                    for u in (int(m.group(1)), int(m.group(2))):
                        out |= support(u, depth + 1)
            elif kind == 'free':
                out |= support(t, depth + 1)
    memo[w] = frozenset(out)
    return memo[w]


cands = set()
for a in MUX.E.res:
    for u in re.findall(r'x(\d+)', a): cands.add(int(u))
for w in sorted(cands):
    try: support(w)
    except Exception: pass

ON_EXP = [3, 10]
ON = set(exp2sel[e] for e in ON_EXP)
v = drive(ON)


def predict(sup):
    live = sorted(sel2exp[s] for s in (sup & ON))
    if not live: return None
    R = FD.INF
    for e in live: R = FD.add(R, bypow[e])
    return R


rows = []
for w, sup in memo.items():
    if w in leafsel or not sup: continue
    live = sup & ON
    if not live: continue
    pr = predict(sup)
    if pr is None: continue
    xv = (v[w] + S) % P
    okx = (xv == pr[0]); oky = (v[w] == pr[1])
    rows.append((len(sup), len(live), w, okx, oky))
rows.sort()
print('ON exponents', ON_EXP, ' wires with live support:', len(rows))
print('%-7s %-6s %-8s %-6s %-6s' % ('|supp|', 'live', 'wire', 'x==pred', 'y==pred'))
shown = 0
for n, nl, w, okx, oky in rows:
    if shown < 30:
        print('%-7d %-6d x%-7d %-6s %-6s' % (n, nl, w, okx, oky))
        shown += 1
agree = [r for r in rows if r[3] or r[4]]
print('\nwires where the model matches on at least one coordinate: %d / %d' % (len(agree), len(rows)))
first_multi = [r for r in rows if r[1] >= 2]
print('wires whose support has >=2 live leaves: %d' % len(first_multi))
if first_multi:
    n, nl, w, okx, oky = first_multi[0]
    print('shallowest multi-live wire: x%d  |supp|=%d live=%d  x_ok=%s y_ok=%s' % (w, n, nl, okx, oky))
    print('  its support exponents:', sorted(sel2exp[s] for s in memo[w]))
    print('  live ones:', sorted(sel2exp[s] for s in (memo[w] & ON)))
    print('  derived by:', C.names[C.trace[w]][:110] if C.trace.get(w) is not None else 'SEEDED')
