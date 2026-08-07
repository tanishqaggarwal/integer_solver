"""Agent B: recover the circuit orientation.

Rule: a variable occurring in exactly ONE gate must be that gate's output.
Then extend to a maximum matching, then break cycles to get a DAG.
"""
import pickle, collections, json, sys, time

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
NV = 38748
NF = len(facs)

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
    cands.append(sorted(allv - sq))
candset = [set(c) for c in cands]

def main():
    matchF = [-1]*NF     # gate -> output var
    matchV = [-1]*NV     # var -> definer gate
    # seed: degree-1 vars
    seeds = 0
    for v in range(NV):
        if len(occ[v]) == 1:
            f = occ[v][0]
            if v in candset[f] and matchF[f] == -1:
                matchF[f] = v; matchV[v] = f; seeds += 1
    print("seeded from degree-1 vars:", seeds)
    # augment with Kuhn's algorithm (iterative DFS)
    def try_aug(start):
        visited = set()
        stack = [(start, iter(cands[start]))]
        path = []
        while stack:
            g, it = stack[-1]
            adv = False
            for v in it:
                if v in visited: continue
                visited.add(v)
                g2 = matchV[v]
                if g2 == -1:
                    # augment along path
                    path.append((g, v))
                    for gg, vv in reversed(path):
                        matchF[gg] = vv; matchV[vv] = gg
                    return True
                stack.append((g2, iter(cands[g2])))
                path.append((g, v))
                adv = True
                break
            if not adv:
                stack.pop()
                if path: path.pop()
        return False
    n = seeds
    t0 = time.time()
    for f in range(NF):
        if matchF[f] == -1 and cands[f]:
            if try_aug(f): n += 1
    print("matching size %d  (%.1fs)" % (n, time.time()-t0))
    print("unmatched gates (assertions):", sum(1 for f in range(NF) if matchF[f] == -1))
    print("unmatched vars (free inputs):", sum(1 for v in range(NV) if matchV[v] == -1))
    # build DAG: output var -> input vars, check acyclicity by topological sort
    indeg = [0]*NF
    users = collections.defaultdict(list)   # var -> gates that USE it (not their output)
    for f in range(NF):
        if matchF[f] == -1: continue
        for v in fvars[f]:
            if v != matchF[f]:
                users[v].append(f)
                indeg[f] += 1
    free = set(v for v in range(NV) if matchV[v] == -1)
    ready = collections.deque()
    known = set(free)
    cnt = [0]*NF
    for f in range(NF):
        if matchF[f] == -1: continue
        cnt[f] = sum(1 for v in fvars[f] if v != matchF[f] and v not in known)
        if cnt[f] == 0: ready.append(f)
    order = []
    while ready:
        f = ready.popleft(); order.append(f)
        v = matchF[f]; known.add(v)
        for g in users[v]:
            cnt[g] -= 1
            if cnt[g] == 0: ready.append(g)
    print("topologically evaluable gates: %d of %d matched" % (len(order), n))
    print("cyclic (unevaluable) matched gates:", n - len(order))
    print("vars reachable:", len(known), "of", NV)
    pickle.dump({'matchF': matchF, 'matchV': matchV, 'order': order, 'free': sorted(free)},
                open(W+'orient.pkl','wb'), -1)
    print("wrote orient.pkl")

if __name__ == '__main__':
    main()
