"""Agent B: post-pass that converts assertion gates into definitions where a single
free input can absorb them, keeping the evaluation order acyclic."""
import pickle, collections, sys, os

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs = M['facs']
NV = 38748; NF = len(facs)
fvars = []
for p in facs:
    vs = set()
    for m in p: vs.update(m)
    fvars.append(vs)
cands = []
for p in facs:
    sq = set()
    for m in p:
        if len(m) == 2 and m[0] == m[1]: sq.add(m[0])
    allv = set()
    for m in p: allv.update(m)
    cands.append(allv - sq)

def main():
    src = sys.argv[1]; dst = sys.argv[2]
    O = pickle.load(open(W+src, 'rb'))
    matchF = list(O['matchF']); matchV = list(O['matchV'])
    order = list(O['order']); free = set(O['free']); assertions = list(O['assertions'])
    excl = set(O['excl'])
    rounds = 0
    while True:
        pos = [-1]*NV
        for k, f in enumerate(order): pos[matchF[f]] = k
        # earliest use of each var in the order
        firstuse = [len(order)+1]*NV
        for k, f in enumerate(order):
            for v in fvars[f]:
                if v != matchF[f] and k < firstuse[v]: firstuse[v] = k
        moved = 0
        for f in list(assertions):
            if f in excl: continue
            U = [v for v in fvars[f] if v in free and v in cands[f]]
            if len(U) != 1: continue
            v = U[0]
            others = [u for u in fvars[f] if u != v]
            if any(u in free for u in others): continue
            ip = max([pos[u] for u in others] + [-1])
            if firstuse[v] <= ip: continue         # would create a cycle
            order.insert(ip+1, f)
            matchF[f] = v; matchV[v] = f
            free.discard(v); assertions.remove(f)
            moved += 1
            # refresh pos/firstuse lazily: restart the sweep
            break
        rounds += 1
        if moved == 0: break
        if rounds % 50 == 0:
            print("  round %d  free=%d assertions=%d" % (rounds, len(free), len(assertions)), flush=True)
    print("post-pass: DEFINED=%d FREE=%d ASSERTIONS=%d" % (len(order), len(free), len(assertions)))
    pickle.dump({'matchF': matchF, 'matchV': matchV, 'order': order,
                 'free': sorted(free), 'assertions': assertions,
                 'live': O.get('live', []), 'excl': sorted(excl)}, open(W+dst,'wb'), -1)
    print("wrote", dst)

if __name__ == '__main__':
    main()
