"""Agent B: symbolic circuit recovery.

Greedily build a topological order: a factor whose unknown-var set is a single
variable v (appearing linearly) DEFINES v.  When stuck, promote the most useful
unknown variable to FREE INPUT.  Result: a DAG + a free-input set + assertions.
"""
import pickle, collections, sys, time, json

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']
NV = 38748

fvars = []
linvars = []   # vars appearing ONLY in linear monomials of the factor
for p in facs:
    s = set(); nlin = set()
    for m in p:
        s.update(m)
        if len(m) > 1: nlin.update(m)
    fvars.append(s); linvars.append(s - nlin)
occ = collections.defaultdict(set)
for i, s in enumerate(fvars):
    for v in s: occ[v].add(i)

def main():
    known = set()
    definer = {}          # var -> factor
    defines = {}          # factor -> var
    assertions = []       # factors that closed with no new var
    free = []
    unresolved = set(range(len(facs)))
    t0 = time.time()
    # count of unknown vars per factor
    nunk = {i: len(fvars[i]) for i in range(len(facs))}
    ready = collections.deque(i for i in unresolved if nunk[i] <= 1)
    rounds = 0
    while unresolved:
        progressed = False
        while ready:
            i = ready.popleft()
            if i not in unresolved: continue
            rem = fvars[i] - known
            if len(rem) > 1: continue
            unresolved.discard(i)
            progressed = True
            if not rem:
                assertions.append(i); continue
            v = next(iter(rem))
            if v in linvars[i]:
                definer[v] = i; defines[i] = v; known.add(v)
                for j in occ[v]:
                    if j in unresolved:
                        nunk[j] -= 1
                        if nunk[j] <= 1: ready.append(j)
            else:
                # v only appears nonlinearly -> cannot define; treat as assertion-ish
                unresolved.add(i)
        if not unresolved: break
        if not progressed:
            # stuck: promote a free input
            cnt = collections.Counter()
            for i in unresolved:
                for v in fvars[i] - known: cnt[v] += 1
            if not cnt: break
            v = cnt.most_common(1)[0][0]
            free.append(v); known.add(v)
            for j in occ[v]:
                if j in unresolved:
                    nunk[j] -= 1
                    if nunk[j] <= 1: ready.append(j)
            rounds += 1
            if rounds % 200 == 0:
                print("  free=%d known=%d unresolved=%d %.1fs" % (len(free), len(known), len(unresolved), time.time()-t0), flush=True)
    print("DONE: defined=%d free=%d assertions=%d unresolved=%d known=%d  %.1fs" %
          (len(defines), len(free), len(assertions), len(unresolved), len(known), time.time()-t0))
    undef = [v for v in range(NV) if v not in known]
    print("vars never known:", len(undef))
    pickle.dump({'definer': definer, 'defines': defines, 'assertions': assertions,
                 'free': free, 'unresolved': list(unresolved)}, open(W+'dag2.pkl','wb'), -1)
    print("wrote dag2.pkl")

if __name__ == '__main__':
    main()
