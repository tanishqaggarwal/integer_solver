#!/usr/bin/env python3
"""One recorded mod-p run: saves val, reason, decision order.

Faster decision heuristic (incremental counter of 2-unknown atoms per var).
"""
import pickle, os, collections, sys, random, time, json
from model import Model, load_assign
from fp import FpEngine, P, sqrt_p

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


def run(policy='wit', seed=1, blacklist=None, verbose=True, tag='run'):
    M = Model(); E = FpEngine(M)
    wit = load_assign(os.path.join(HERE, '..', 'best',
                                   'new_instance_partial_39026.json'))
    rng = random.Random(seed)
    bl = set(blacklist or [])
    val = [None] * NV
    reason = [None] * NV
    order = []
    unk = [len(E.avarlist[a]) for a in range(M.na)]
    two = collections.Counter()          # var -> #atoms with exactly 2 unknowns
    for a in range(M.na):
        if unk[a] == 2:
            for x in E.avarlist[a]:
                two[x] += 1
    Q = collections.deque(a for a in range(M.na) if unk[a] <= 1)
    decisions = []
    t0 = time.time()

    def assign(u, x, r):
        val[u] = x; reason[u] = r; order.append(u)
        for b in E.var2atoms[u]:
            before = unk[b]
            unk[b] = before - 1
            if before == 2:
                for y in E.avarlist[b]:
                    if val[y] is None:
                        two[y] -= 1
            elif before == 3:
                for y in E.avarlist[b]:
                    if val[y] is None and y != u:
                        two[y] += 1
            if unk[b] <= 1:
                Q.append(b)

    def propagate():
        branch = []
        while Q:
            a = Q.popleft()
            miss = [x for x in E.avarlist[a] if val[x] is None]
            if len(miss) > 1:
                continue
            if not miss:
                if E.eval_atom(a, val) != 0:
                    return a, branch
                continue
            u = miss[0]
            c2, c1, c0 = E.reduce(a, val, u)
            if c2 == 0:
                if c1 == 0:
                    if c0 != 0:
                        return a, branch
                    continue
                assign(u, (-c0) * pow(c1, -1, P) % P, a)
            else:
                d = (c1 * c1 - 4 * c2 * c0) % P
                s = sqrt_p(d)
                if s is None:
                    return a, branch
                inv = pow(2 * c2 % P, -1, P)
                roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
                if len(roots) == 1:
                    assign(u, roots[0], a)
                else:
                    branch.append((u, a, roots))
        return None, branch

    def cone_of(a):
        seen = set(); decs = set()
        st = list(E.avarlist[a])
        while st:
            v = st.pop()
            if v in seen or val[v] is None:
                continue
            seen.add(v)
            r = reason[v]
            if r is None or r == 'dec':
                decs.add(v)
            else:
                st.extend(E.avarlist[r])
        return decs, seen

    status = 'ok'; badatom = None; decs = None; cone = None
    while True:
        bad, branch = propagate()
        if bad is not None:
            decs, cone = cone_of(bad)
            status = 'conflict'; badatom = bad
            break
        pend = [(u, a, r) for u, a, r in branch if val[u] is None]
        if pend:
            for u, a, roots in pend:
                if val[u] is not None:
                    continue
                if policy == 'wit':
                    w = wit[u] % P; x = w if w in roots else roots[0]
                elif policy == 'zero':
                    x = 0 if 0 in roots else roots[0]
                else:
                    x = rng.choice(roots)
                assign(u, x, 'dec'); decisions.append(u)
            continue
        cands = [(n, v) for v, n in two.items() if n > 0 and val[v] is None and v not in bl]
        if cands:
            u = max(cands)[1]
        else:
            rest = [v for v in range(NV) if val[v] is None]
            if not rest:
                break
            r2 = [v for v in rest if v not in bl]
            u = r2[0] if r2 else rest[0]
        if policy == 'wit':
            x = wit[u] % P
        elif policy == 'zero':
            x = 0
        else:
            x = rng.randrange(P)
        assign(u, x, 'dec'); decisions.append(u)
        if verbose and len(decisions) % 500 == 0:
            print(f"  dec={len(decisions)} known={len(order)} t={time.time()-t0:.0f}s", flush=True)
    known = sum(1 for x in val if x is not None)
    print(f"[{tag}] status={status} bad={badatom} known={known}/{NV} "
          f"decisions={len(decisions)} t={time.time()-t0:.0f}s", flush=True)
    if status == 'conflict':
        print(f"   atom: {M.src[badatom][:140]}")
        print(f"   cone={len(cone)} decisions_in_cone={len(decs)}: {sorted(decs)[:30]}")
    out = {'val': val, 'reason': reason, 'order': order,
           'decisions': decisions, 'status': status, 'bad': badatom,
           'decs': sorted(decs) if decs else None}
    pickle.dump(out, open(os.path.join(HERE, f'fprun_{tag}.pkl'), 'wb'))
    return out, M, E


if __name__ == '__main__':
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    sd = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    blp = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, 'blacklist.json')
    bl = json.load(open(blp)) if os.path.exists(blp) else []
    run(pol, sd, bl, tag=f'{pol}_{sd}')
