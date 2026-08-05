"""Branch-and-bound: is there ANY relaxation set with failing-equation cost <= 10
that makes the chain-1-relaxed mod-p system consistent?  (cost 9 => 39,024, 10 => 39,023)"""
import pickle, collections, sys, os, time, heapq
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0,S9); os.chdir(S9)
P = 2**256-2**32-977
d = pickle.load(open('chain2/jac24.pkl','rb')); J = d['J']; base = d['base']
a2e = {a:set(v) for a,v in pickle.load(open('atom2eq.pkl','rb')).items()}
BASECOST = a2e.get(22229, set())          # chain 1 relaxed
LIMIT = int(sys.argv[1]) if len(sys.argv)>1 else 10

ROWS_ALL = collections.defaultdict(dict)
for f, col in J.items():
    for a, dv in col.items():
        x = dv % P
        if x: ROWS_ALL[a][f] = x
RHS_ALL = {a: (-base.get(a,0)) % P for a in ROWS_ALL}

def eliminate(drop):
    rows = {a: r for a, r in ROWS_ALL.items() if a not in drop}
    order = sorted(rows, key=lambda a: len(rows[a]))
    pivots = {}
    for a in order:
        row = dict(rows[a]); r = RHS_ALL[a]; pv = {a: 1}
        while True:
            hit = next((c for c in row if c in pivots), None)
            if hit is None: break
            prow, prhs, ppv = pivots[hit]
            fac = row[hit]*pow(prow[hit], P-2, P) % P
            for c, val in prow.items():
                nv = (row.get(c,0) - fac*val) % P
                if nv: row[c] = nv
                elif c in row: del row[c]
            r = (r - fac*prhs) % P
            for k, val in ppv.items(): pv[k] = (pv.get(k,0) - fac*val) % P
        if not row:
            if r: return {k: c for k, c in pv.items() if c}
            continue
        c = min(row); pivots[c] = (row, r, pv)
    return None

seen = set(); t0 = time.time(); nodes = 0
pq = [(len(BASECOST), tuple(), BASECOST)]
best = None
while pq:
    cost, drop, union = heapq.heappop(pq)
    key = frozenset(drop)
    if key in seen: continue
    seen.add(key); nodes += 1
    cert = eliminate(set(drop))
    if cert is None:
        best = (cost, drop); print(f'\n*** CONSISTENT at cost {cost} with drop={sorted(drop)}')
        break
    cands = []
    for a in cert:
        if a in drop: continue
        nu = union | a2e.get(a, set())
        if len(nu) <= LIMIT: cands.append((len(nu), a, nu))
    cands.sort()
    for nc, a, nu in cands[:14]:
        nd = tuple(sorted(drop + (a,)))
        if frozenset(nd) not in seen: heapq.heappush(pq, (nc, nd, nu))
    if nodes % 25 == 0:
        print(f'  nodes={nodes} frontier={len(pq)} best_cost_seen={cost} t={time.time()-t0:.0f}s')
    if time.time()-t0 > 1500: print('time limit'); break
print(f'\nnodes explored: {nodes}, time {time.time()-t0:.0f}s')
if best is None:
    print(f'NO relaxation set with cost <= {LIMIT} makes the system consistent.')
else:
    pickle.dump(best, open('chain2/bb_best.pkl','wb'))
