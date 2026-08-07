#!/usr/bin/env python3
"""No-good loop over the *choice of free inputs* in the mod-p solve.

Each round: propagate; if stuck, designate a free input (avoiding blacklisted
variables); on conflict, blacklist the non-boolean decision variables in the
conflict cone and restart.  The blacklist is persisted to disk after every
round so a resumed run keeps the pruning it paid for.
"""
import pickle, os, collections, sys, random, time, json
from model import Model, load_assign
from fp import FpEngine, P, sqrt_p

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748
BLPATH = os.path.join(HERE, 'blacklist.json')


def load_bl():
    if os.path.exists(BLPATH):
        return set(json.load(open(BLPATH)))
    return set()


def save_bl(bl):
    json.dump(sorted(bl), open(BLPATH, 'w'))


def one_run(M, E, wit, policy, rng, blacklist, verbose=False):
    val = [None] * NV
    reason = [None] * NV
    unk = [len(E.avarlist[a]) for a in range(M.na)]
    Q = collections.deque(a for a in range(M.na) if unk[a] <= 1)
    decisions = []

    def assign(u, x, r):
        val[u] = x; reason[u] = r
        for b in E.var2atoms[u]:
            unk[b] -= 1
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

    while True:
        bad, branch = propagate()
        if bad is not None:
            decs, cone = cone_of(bad)
            return ('conflict', bad, decs, cone, val, reason, decisions)
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
        unknown = [v for v in range(NV) if val[v] is None]
        if not unknown:
            return ('ok', None, None, None, val, reason, decisions)
        deg = collections.Counter()
        for a in range(M.na):
            if unk[a] == 2:
                for x in E.avarlist[a]:
                    if val[x] is None:
                        deg[x] += 1
        cands = [(n, v) for v, n in deg.items() if v not in blacklist]
        if cands:
            u = max(cands)[1]
        else:
            free_ok = [v for v in unknown if v not in blacklist]
            u = free_ok[0] if free_ok else unknown[0]
        if policy == 'wit':
            x = wit[u] % P
        elif policy == 'zero':
            x = 0
        else:
            x = rng.randrange(P)
        assign(u, x, 'dec'); decisions.append(u)


def main():
    policy = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    M = Model(); E = FpEngine(M)
    wit = load_assign(os.path.join(HERE, '..', 'best',
                                   'new_instance_partial_39026.json'))
    rng = random.Random(seed)
    bl = load_bl()
    boolvars = set()
    for a in range(M.na):
        q = M.polys[a]
        vs = set()
        for m in q:
            vs |= set(m)
        if len(vs) == 1 and max(len(m) for m in q) == 2:
            boolvars |= vs
    print("boolean-constrained vars:", len(boolvars), "blacklist start:", len(bl))
    t0 = time.time()
    for it in range(rounds):
        res = one_run(M, E, wit, policy, rng, bl)
        kind = res[0]
        if kind == 'ok':
            val = res[4]
            print(f"[{it}] COMPLETE mod-p solution! decisions={len(res[6])} t={time.time()-t0:.0f}s")
            pickle.dump(val, open(os.path.join(HERE, f'fp_full_{policy}_{seed}.pkl'), 'wb'))
            return
        _, bad, decs, cone, val, reason, decisions = res
        new = sorted(v for v in decs if v not in boolvars and v not in bl)
        known = sum(1 for x in val if x is not None)
        print(f"[{it}] conflict a{bad} ({M.src[bad][:70]}) cone={len(cone)} "
              f"decs={len(decs)} new_bl={len(new)} known={known} t={time.time()-t0:.0f}s",
              flush=True)
        if not new:
            print("   no new blacklist candidates -> stuck.  decs:", sorted(decs)[:20])
            print("   boolean decisions in cone:",
                  sorted(v for v in decs if v in boolvars)[:20])
            break
        bl |= set(new)
        save_bl(bl)
    save_bl(bl)


if __name__ == '__main__':
    main()
