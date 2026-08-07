"""Agent B: ACYCLIC circuit orientation by topological growth.
Free inputs are promoted only when the frontier stalls (chosen by max unblocking)."""
import pickle, collections, heapq, time, sys

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

def main():
    excl = set()
    if len(sys.argv) > 1:
        excl = set(pickle.load(open(sys.argv[1],'rb')))
        print("excluding %d gates from being definitions" % len(excl))
    known = set()
    matchF = [-1]*NF
    matchV = [-1]*NV
    order = []
    assertions = []
    free = []
    nunk = [len(fvars[i]) for i in range(NF)]
    live = [True]*NF
    for i in excl: live[i] = False
    ready = collections.deque(i for i in range(NF) if nunk[i] <= 1 and live[i])
    # lazy max-heap on how many live gates a var blocks
    blk = collections.Counter()
    for i in range(NF):
        if not live[i]: continue
        for v in fvars[i]: blk[v] += 1
    heap = [(-n, v) for v, n in blk.items()]
    heapq.heapify(heap)
    t0 = time.time()
    def learn(v):
        known.add(v)
        for g in occ[v]:
            if live[g]:
                nunk[g] -= 1
                blk[v] -= 1
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
            if v not in cands[i]:
                continue   # cannot solve for it here; leave the gate live
            live[i] = False
            matchF[i] = v; matchV[v] = i; order.append(i)
            for u in fvars[i]: blk[u] -= 1
            learn(v)
        # stall
        rest = [i for i in range(NF) if live[i]]
        if not rest: break
        # pick the unknown var blocking the most live gates
        best = None
        while heap:
            neg, v = heap[0]
            if v in known:
                heapq.heappop(heap); continue
            cur = sum(1 for g in occ[v] if live[g])
            if -neg > cur:
                heapq.heapreplace(heap, (-cur, v)); continue
            best = v; break
        if best is None: break
        free.append(best)
        heapq.heappop(heap)
        learn(best)
        if len(free) % 500 == 0:
            print("  free=%d known=%d live=%d %.1fs" % (len(free), len(known), sum(live), time.time()-t0), flush=True)
    print("DEFINED=%d  FREE=%d  ASSERTIONS=%d  live-left=%d  known=%d  %.1fs" %
          (len(order), len(free), len(assertions), sum(live), len(known), time.time()-t0))
    pickle.dump({'matchF': matchF, 'matchV': matchV, 'order': order,
                 'free': free, 'assertions': assertions,
                 'live': [i for i in range(NF) if live[i]], 'excl': sorted(excl)}, open(W+(sys.argv[2] if len(sys.argv)>2 else 'orient2.pkl'),'wb'), -1)
    print("wrote orient2.pkl")

if __name__ == '__main__':
    main()
