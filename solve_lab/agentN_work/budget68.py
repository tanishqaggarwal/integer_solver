"""Collateral budget |W| = 1 on the COMPLETE 68-knob set, exact polynomial model, resumable.

To beat 39,026 at budget W the region gain must satisfy (12 - g) + W < 7, i.e. g >= W + 6.
At W = 0 that is g >= 6, refuted and mod-p certified.  Here W = 1: every one of the collateral rows
is dropped in turn, the saturation loop is redone exactly on the enlarged constraint set, and the
maximum integrally zeroable region-row count is computed with zsolve.  Writes one JSON line per
drop so it can be resumed after a restart.
"""
import os, sys, json, time, pickle
from collections import defaultdict
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN, zsolve
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits
from polyexact import P
from polyfull import exact_polys
from kerquad import int_kernel_columns

OUT = os.path.join(HERE, 'runs', 'budget68.jsonl')


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
print('68 knobs, %d live collateral rows, region %d' % (len(live), len(Rl)), flush=True)


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


def reduce_and_solve(dropset):
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
            return 0, 0, [], 'lattice collapsed'
        Kn = lll(Kn)
        K = lll([[sum(u[a] * K[a][j] for a in range(len(u))) for j in range(k)] for u in Kn])
        live_o = nl
    # remaining nonlinear collateral rows: only pure squares are handled exactly here
    dead = set()
    other = []
    for e in live_o:
        c = restrict(polys[e], K).c
        if len(c) == 1:
            mono = list(c)[0]
            idx = [j for j, x in enumerate(mono) if x]
            if len(idx) == 1 and mono[idx[0]] == 2:
                dead.add(idx[0])
                continue
        if c:
            other.append(e)
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
    if other:
        return -2, len(freec), [], 'unhandled nonlinear rows %s' % other
    opt, rws, exh, nd = zsolve.max_zero_rows(M, b, len(freec), len(M), node_cap=2000000)
    return opt, len(freec), [Rl[i] for i in rws], 'ok' if exh else 'capped'


done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        try:
            done.add(json.loads(line)['drop'])
        except Exception:
            pass
print('already done: %d' % len(done), flush=True)

t00 = time.time()
best = (5, None)
with open(OUT, 'a') as f:
    for i, e in enumerate(live):
        if e in done:
            continue
        t0 = time.time()
        opt, rk, rws, note = reduce_and_solve({e})
        rec = dict(drop=e, opt=opt, rank=rk, rows=rws, note=note,
                   failing=(12 - opt + 1) if opt >= 0 else None, secs=round(time.time() - t0, 1))
        f.write(json.dumps(rec) + '\n')
        f.flush()
        if opt > best[0]:
            best = (opt, e)
            print('   NEW BEST g=%d dropping eq %d (score %d)'
                  % (opt, e, 39033 - (12 - opt + 1)), flush=True)
        if i % 20 == 0:
            print('   %d/%d  (%.0fs elapsed)  best g so far %d'
                  % (i, len(live), time.time() - t00, best[0]), flush=True)
print('DONE |W|=1 on 68 knobs: best g = %d (need >= 7 to beat 39,026)' % best[0], flush=True)
