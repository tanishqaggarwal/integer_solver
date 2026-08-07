#!/usr/bin/env python3
"""K23: find the leaf the mux descent missed and settle its side by walking the wiring
up from that leaf, not by guessing."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
import mux as MUX
from parse import node_str

rs = json.load(open(K + '/rootsupport.json'))
ch = json.load(open(K + '/chain.json'))
sel2exp = {}
exp2sel = {}
for i in range(256):
    s = ch['sel'][str(i)]; e = ch['exp'][str(i)]
    sel2exp[s] = e; exp2sel[e] = s
IA = set(sel2exp[s] for s in rs['A.x']) | set(sel2exp[s] for s in rs['A.y'])
IB = set(sel2exp[s] for s in rs['B.x']) | set(sel2exp[s] for s in rs['B.y'])
miss = sorted(set(range(256)) - IA - IB)
print('missing exponents:', miss, '-> selectors', [exp2sel[e] for e in miss])

D = FD.points()
wof = {}
for l in D['leaves']:
    wof[l['sel']] = (l['wx'], l['wy'])

# build value-wire -> parent-slot edges over the whole circuit
gatedpat = re.compile(r'^\(x(\d+)\*x(\d+)\)$')
E = MUX.E
allwires = set()
for a in E.res: pass
parent = collections.defaultdict(set)
slots = []
for a in E.res:
    m = re.match(r'^\(\(x(\d+)-x(\d+)\)-', a) or re.match(r'^\(\((\d+)\*\(x(\d+)-x(\d+)\)\)-', a)
    pass
# simpler: for every variable w that has a decodable source, record its branch children
cands = set()
for a in E.res:
    for u in re.findall(r'x(\d+)', a): cands.add(int(u))
print('candidate wires:', len(cands))
edges = 0
for w in sorted(cands):
    try:
        src = MUX.source_of(w)
    except Exception:
        continue
    for z, coef in src:
        if z == 'CONST': continue
        for kind, t in MUX.mux_terms(z):
            if kind == 'gated':
                mm = gatedpat.match(t)
                if mm:
                    for u in (int(mm.group(1)), int(mm.group(2))):
                        parent[u].add(w); edges += 1
            elif kind == 'free':
                parent[t].add(w); edges += 1
print('parent edges:', edges)

TARGETS = {12186: 'A.x', 16742: 'A.y', 14853: 'B.x', 24908: 'B.y'}
for e in miss:
    s = exp2sel[e]
    for w in wof[s]:
        seen = set([w]); st = [w]; hit = set()
        while st:
            u = st.pop()
            if u in TARGETS: hit.add(TARGETS[u]); continue
            for q in parent.get(u, ()):
                if q not in seen: seen.add(q); st.append(q)
        print('exponent', e, 'selector x%d wire x%d -> reaches' % (s, w), sorted(hit), 'visited', len(seen))
