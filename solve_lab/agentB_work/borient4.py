"""Agent B: rule-based + topological circuit orientation (v4).

Rule for a gate containing a degree-2 monomial:
  * if some variable inside a quadratic monomial occurs in exactly ONE gate, it is a
    free QUOTIENT HANDLE and is the output;
  * else the output is the variable in linear position (the product's result).
Pure linear gates keep the greedy topological choice.
"""
import pickle, collections, heapq, time, sys, json

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
NV = 38748; NF = len(facs)

occ = collections.defaultdict(list)
fvars = []
for i, p in enumerate(facs):
    vs = set()
    for m in p: vs.update(m)
    fvars.append(vs)
    for v in vs: occ[v].append(i)
cands = []
for p in facs:
    sq = set()
    for m in p:
        if len(m) == 2 and m[0] == m[1]: sq.add(m[0])
    allv = set()
    for m in p: allv.update(m)
    cands.append(allv - sq)

pref = [None]*NF
for i, p in enumerate(facs):
    quad = set(); lin = set()
    for m, c in p.items():
        if len(m) > 1: quad.update(m)
        elif len(m) == 1: lin.add(m[0])
    lin -= quad
    if not quad: continue
    h = [v for v in quad if len(occ[v]) == 1 and v in cands[i]]
    if len(h) == 1: pref[i] = h[0]; continue
    if len(lin) == 1:
        v = next(iter(lin))
        if v in cands[i]: pref[i] = v
print("gates with a forced output:", sum(1 for x in pref if x is not None))

def main():
    excl = set()
    if len(sys.argv) > 1 and sys.argv[1] != '-':
        excl = set(pickle.load(open(sys.argv[1],'rb')))
        print("excluding %d gates" % len(excl))
    out = sys.argv[2] if len(sys.argv) > 2 else 'orient4.pkl'
    known = set(); matchF = [-1]*NF; matchV = [-1]*NV
    order = []; assertions = []; free = []
    live = [True]*NF
    for i in excl: live[i] = False
    nunk = [len(fvars[i]) for i in range(NF)]
    ready = collections.deque(i for i in range(NF) if live[i] and nunk[i] <= 1)
    blk = collections.Counter()
    for i in range(NF):
        if live[i]:
            for v in fvars[i]: blk[v] += 1
    heap = [(-n, v) for v, n in blk.items()]; heapq.heapify(heap)
    t0 = time.time()
    def learn(v):
        known.add(v)
        for g in occ[v]:
            if live[g]:
                nunk[g] -= 1
                if nunk[g] <= 1: ready.append(g)
    while True:
        while ready:
            i = ready.popleft()
            if not live[i]: continue
            rem = fvars[i] - known
            if len(rem) > 1: continue
            if not rem:
                live[i] = False; assertions.append(i)
                for v in fvars[i]: blk[v] -= 1
                continue
            v = next(iter(rem))
            if v not in cands[i]: continue
            if pref[i] is not None and pref[i] != v:
                continue          # wrong direction; leave the gate as a constraint
            live[i] = False
            matchF[i] = v; matchV[v] = i; order.append(i)
            for u in fvars[i]: blk[u] -= 1
            learn(v)
        rest = [i for i in range(NF) if live[i]]
        if not rest: break
        best = None
        while heap:
            neg, v = heap[0]
            if v in known: heapq.heappop(heap); continue
            cur = sum(1 for g in occ[v] if live[g])
            if -neg > cur: heapq.heapreplace(heap, (-cur, v)); continue
            best = v; break
        if best is None:
            # nothing left to promote: remaining live gates that are fully known are assertions
            for i in rest:
                if not (fvars[i] - known):
                    live[i] = False; assertions.append(i)
            if not any(live): break
            # promote any unknown var
            un = set()
            for i in range(NF):
                if live[i]: un |= (fvars[i] - known)
            if not un:
                for i in range(NF):
                    if live[i]: live[i] = False; assertions.append(i)
                break
            best = min(un)
        else:
            heapq.heappop(heap)
        free.append(best); learn(best)
    print("DEFINED=%d FREE=%d ASSERTIONS=%d live=%d known=%d %.1fs" %
          (len(order), len(free), len(assertions), sum(live), len(known), time.time()-t0))
    pickle.dump({'matchF': matchF, 'matchV': matchV, 'order': order,
                 'free': free, 'assertions': assertions,
                 'live': [i for i in range(NF) if live[i]], 'excl': sorted(excl)},
                open(W+out,'wb'), -1)
    print("wrote", out)

if __name__ == '__main__':
    main()
