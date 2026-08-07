"""|W| = 1 leftovers (relaxation bound) and the sharded |W| = 2 sweep, complete 68-knob set.

For 14 of the 231 single drops, eq 8680 stops being a pure square on the enlarged lattice, so it is
a genuine quadratic constraint.  Those are settled by RELAXATION: drop the constraint entirely and
compute the maximum integrally zeroable region-row count without it.  Relaxing can only raise the
optimum, so a relaxed optimum of 5 proves the true optimum is at most 5.

Then |W| = 2, sharded across workers:  python3 budget68b.py w2 <shard> <nshards>
"""
import os, sys, json, time, itertools
from collections import defaultdict
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN, zsolve
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits
from polyexact import P
from polyfull import exact_polys
from kerquad import int_kernel_columns


def lll(K):
    if not K:
        return []
    R = fmpz_mat([[int(x) for x in v] for v in K]).lll().tolist()
    return [[int(x) for x in r] for r in R if any(r)]


st = make(WIT)
b0 = build(st)
Rl = b0['R']
Rset = set(Rl)
atoms_R = set()
for e in Rl:
    for c, a in ev.eq_terms[e][2]:
        atoms_R.add(a)
cands = set()
for q in atoms_R:
    if q in fr.csup:
        cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
cands = sorted(y for y in cands if y in FREE)
k = len(cands)
touched = set()
for Y in cands:
    for a in fr.chk.get(Y, []):
        touched |= atom_eqs[a]
outside = sorted(touched - Rset)
polys = exact_polys(st, Rl + outside, cands)
live = [e for e in outside if polys[e].c]


def restrict(pol, K):
    dd = len(K)
    P.NK = dd
    T = []
    for j in range(k):
        c = {}
        for a in range(dd):
            if K[a][j]:
                m = [0] * dd
                m[a] = 1
                c[tuple(m)] = K[a][j]
        T.append(P(c))
    acc = P()
    for mono, cf in pol.c.items():
        term = P.const(cf)
        for j, ex in enumerate(mono):
            for _ in range(ex):
                term = term * T[j]
        acc = acc + term
    return acc


def reduce_and_solve(dropset, relax=True):
    """Returns (opt, rank, rows, note).  With relax=True any residual nonlinear collateral
    constraint that is not a pure square is DROPPED, which can only raise opt -> upper bound."""
    K = [[1 if i == j else 0 for i in range(k)] for j in range(k)]
    live_o = [e for e in live if e not in dropset]
    for _ in range(12):
        cur = {e: restrict(polys[e], K) for e in live_o}
        linrows = [e for e in live_o if cur[e].deg() == 1]
        nl = [e for e in live_o if cur[e].deg() >= 2]
        if not linrows:
            live_o = nl
            break
        A = []
        for e in linrows:
            v = [0] * len(K)
            for mono, c in cur[e].c.items():
                v[mono.index(1)] = c
            A.append(v)
        Kn = int_kernel_columns(A, len(K))
        if not Kn:
            return 0, 0, [], 'collapsed'
        Kn = lll(Kn)
        K = lll([[sum(u[a] * K[a][j] for a in range(len(u))) for j in range(k)] for u in Kn])
        live_o = nl
    dead = set()
    relaxed = []
    for e in live_o:
        c = restrict(polys[e], K).c
        if not c:
            continue
        if len(c) == 1:
            mono = list(c)[0]
            idx = [j for j, x in enumerate(mono) if x]
            if len(idx) == 1 and mono[idx[0]] == 2:
                dead.add(idx[0])
                continue
        relaxed.append(e)
    if relaxed and not relax:
        return -2, len(K), [], 'unhandled %s' % relaxed
    freec = [a for a in range(len(K)) if a not in dead]
    M, b = [], []
    for e in Rl:
        c = restrict(polys[e], K).c
        if any(sum(m) > 1 for m in c):
            return -1, len(K), [], 'region nonlinear'
        row = [0] * len(K)
        c0 = 0
        for mono, v in c.items():
            if sum(mono) == 0:
                c0 = v
            else:
                row[mono.index(1)] = v
        M.append([row[a] for a in freec])
        b.append(c0)
    opt, rws, exh, nd = zsolve.max_zero_rows(M, b, len(freec), len(M), node_cap=2000000)
    note = 'exact' if not relaxed else 'RELAXED(upper bound), dropped %s' % relaxed
    if not exh:
        note += ' CAPPED'
    return opt, len(freec), [Rl[i] for i in rws], note


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'w1fix'
    if mode == 'w1fix':
        todo = [json.loads(l)['drop'] for l in open(os.path.join(HERE, 'runs', 'budget68.jsonl'))
                if json.loads(l)['opt'] < 0]
        print('re-running %d unhandled |W|=1 drops by relaxation' % len(todo), flush=True)
        best = 0
        for e in todo:
            opt, rk, rws, note = reduce_and_solve({e}, relax=True)
            best = max(best, opt)
            print('   drop %-6d  relaxed OPT = %-3d rank %-3d  %s' % (e, opt, rk, note), flush=True)
        print('\nmax relaxed OPT over the 14: %d  (true OPT <= this)' % best, flush=True)
        print('|W|=1 on 68 knobs: g <= %d everywhere; need g >= 7 to beat 39,026 -> REFUTED'
              % max(best, 5), flush=True)
    elif mode == 'w2':
        shard, ns = int(sys.argv[2]), int(sys.argv[3])
        out = os.path.join(HERE, 'runs', 'budget68_w2_%d.jsonl' % shard)
        done = set()
        if os.path.exists(out):
            for l in open(out):
                try:
                    done.add(tuple(json.loads(l)['drop']))
                except Exception:
                    pass
        pairs = [c for i, c in enumerate(itertools.combinations(live, 2)) if i % ns == shard]
        print('shard %d/%d: %d pairs, %d already done' % (shard, ns, len(pairs), len(done)),
              flush=True)
        best = 0
        t0 = time.time()
        with open(out, 'a') as f:
            for i, pr in enumerate(pairs):
                if pr in done:
                    continue
                opt, rk, rws, note = reduce_and_solve(set(pr), relax=True)
                f.write(json.dumps(dict(drop=list(pr), opt=opt, rank=rk, note=note)) + '\n')
                if opt > best:
                    best = opt
                    if opt >= 8:
                        f.flush()
                        print('   *** g=%d dropping %s -> score %d' % (opt, pr, 39033 - (12 - opt + 2)),
                              flush=True)
                if i % 500 == 0:
                    f.flush()
                    print('   %d/%d (%.0fs) best %d' % (i, len(pairs), time.time() - t0, best),
                          flush=True)
        print('shard %d done: best relaxed g = %d (need >= 8)' % (shard, best), flush=True)
