#!/usr/bin/env python3
"""K43: forward-only closure at EVERY slot, then the full validation table.

Every slot wire is pinned by one residual atom.  Any OTHER atom mentioning that wire is a
consumer, and letting the closure solve one of those for the wire is a backward derivation --
the bug that produced two wrong headline results in this directory.  Here every slot wire is
restricted to its own pin, everywhere, not just at the root."""
import sys, os, json, re, collections, time
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from k26_drive import FORBID

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
h, ls, ob, wr = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
dv = [u for u in range(NV) if u not in set(C.E.free)]
ORDER = h + ls + ob + dv + wr
S = FD.SHIFT
ch = json.load(open(K + '/chain.json'))
D = FD.points()
bp = {}
for i_s, e in ch['exp'].items():
    bp[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
e2s = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}

# ---- build the pin map: slot wire -> its own residual pin atom -----------------------
# shapes seen in this file, W is always the slot and Z the mux source:
#   ((xW-xZ)-xH)   ((xW-xZ)+xH)   ((K*(xW-xZ))-xH)   ((xW-xZ)-(K*xH))   ((K*(xW-xZ))-(K*xH))
PAT = [re.compile(r'^\(\(x(\d+)-x(\d+)\)[-+]\(?(?:\d+\*)?x\d+\)?\)$'),
       re.compile(r'^\(\(\d+\*\(x(\d+)-x(\d+)\)\)[-+]\(?(?:\d+\*)?x\d+\)?\)$')]
pin = {}
dup = 0
for i, nm in enumerate(C.names):
    for p in PAT:
        m = p.match(nm)
        if m:
            w = int(m.group(1))
            if w in pin: dup += 1
            else: pin[w] = i
            break
# leaf wires are pinned by their leaf pin, not a slot pin - never let a consumer set them
leafpin = re.compile(r'^\(\(x(\d+)\*\(x(\d+)-\d{20,}\)\)[-+]')
for i, nm in enumerate(C.names):
    m = leafpin.match(nm)
    if m:
        w = int(m.group(2))
        if w not in pin: pin[w] = i
print('slot/leaf wires given an exclusive pin: %d   (ambiguous, skipped: %d)' % (len(pin), dup))


def close(on, guard=True):
    seed = {u: 0 for u in h}
    for u in ls: seed[u] = 1 if u in on else 0
    for u in ob: seed[u] = 0
    v, _ = C.close(seed, ORDER, forbid=FORBID, pin=pin if guard else None)
    return v


def comp(es):
    R = FD.INF
    for e in es: R = FD.add(R, bp[e])
    return R


rs = json.load(open(K + '/rootsupport.json'))
s2e = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
IA = set(s2e[s] for s in rs['A.x']) | set(s2e[s] for s in rs['A.y']) | {163}

TABLE = [[0, 1], [3, 10], [3, 5], [0, 1, 2, 4], [0, 1, 3], [0, 1, 2, 3, 5],
         [4, 6, 7, 10, 12], [0, 2, 4, 6, 8, 3, 5, 10, 13], [5, 10], [3, 5, 10],
         [0], [3], [1, 2], [12, 13]]
print('\n%-28s %-6s %-8s %-8s %-8s %-8s' % ('ON exponents', 'folds', 'A(guard)', 'B(guard)',
                                            'A(no g)', 'B(no g)'))
res = []
for ON in TABLE:
    ea = sorted(e for e in ON if e in IA); eb = sorted(e for e in ON if e not in IA)
    on = set(e2s[e] for e in ON)
    row = []
    for guard in (True, False):
        v = close(on, guard)
        A = ((v[12186] + S) % P, v[16742]); B = ((v[14853] + S) % P, v[24908])
        row += [(A == comp(ea)) if ea else '-', (B == comp(eb)) if eb else '-']
    print('%-28s %-6s %-8s %-8s %-8s %-8s' % (ON, bool(ea) and bool(eb), row[0], row[1], row[2], row[3]))
    res.append((ON, row))

ok = sum(1 for _, r in res for x in r[:2] if x is True)
tot = sum(1 for _, r in res for x in r[:2] if x != '-')
okn = sum(1 for _, r in res for x in r[2:] if x is True)
print('\nGUARDED   : %d / %d halves match' % (ok, tot))
print('UNGUARDED : %d / %d halves match' % (okn, tot))
