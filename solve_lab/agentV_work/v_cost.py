#!/usr/bin/env python3
"""agent V -- V3.  Where the joint solve stops being bounded.

The joint solve enumerates the shift on k-1 wires per prime power and root-finds the last, so its
cost is  SUM over prime powers q^e in the component's moduli  of  q^(e(k-1)).  It is therefore
dominated by the LARGEST prime power Q in the component, and the affordable k is

        k_max(Q, budget)  =  1 + floor( log(budget) / log(Q) )

which is a function of the moduli, NOT of the number of conditions.  This prices every component
of the shares-a-condition graph against a stated budget and reports the distribution, so that the
"there is a size beyond which this stops being bounded" question has a number attached.

Usage: python3 v_cost.py [budget]
"""
import os, sys, json, math, collections
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
SL, p = B.SL, B.p


def comps_structural():
    A2W = {a: sorted(B.wires_of(a)) for a in B.CGT2}
    par = {a: a for a in A2W}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    W2A = collections.defaultdict(list)
    for a, ws in A2W.items():
        for w in ws:
            W2A[w].append(a)
    for w, ats in W2A.items():
        for b in ats[1:]:
            rx, ry = find(ats[0]), find(b)
            if rx != ry:
                par[rx] = ry
    out = collections.defaultdict(list)
    for a in A2W:
        out[find(a)].append(a)
    return [(ats, sorted(set(w for a in ats for w in A2W[a]))) for ats in out.values()]


def biggest_pp(ats):
    Q = 1
    for a in ats:
        for q, e in B.factor(abs(SL[a])//p).items():
            Q = max(Q, q**e)
    return Q


def kmax(Q, budget):
    if Q <= 1:
        return 99
    return 1 + int(math.floor(math.log(budget)/math.log(Q)))


if __name__ == '__main__':
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 10**7
    print('budget = %.0e enumeration points per prime power\n' % budget, flush=True)
    C = comps_structural()
    print('%d structural components over all 927 c>1 conditions\n' % len(C), flush=True)

    rows = []
    for ats, ws in C:
        Q = biggest_pp(ats)
        rows.append((len(ats), len(ws), Q, kmax(Q, budget)))

    print('%-6s %-6s %-12s %-8s %s' % ('#cond', '#wires', 'largest q^e', 'k_max', 'verdict'), flush=True)
    hist = collections.Counter()
    reach = 0
    for nc, nw, Q, km in sorted(rows, key=lambda t: -t[1])[:20]:
        ok = 'FULL COMPONENT REACHABLE' if km >= nw else 'k_max < #wires: only subsets of size %d' % km
        print('%-6d %-6d %-12d %-8d %s' % (nc, nw, Q, km, ok), flush=True)
    for nc, nw, Q, km in rows:
        hist[km] += 1
        if km >= nw:
            reach += 1
    print('\nk_max histogram over the %d components: %s' % (len(rows), dict(sorted(hist.items()))),
          flush=True)
    print('components whose ENTIRE wire set is jointly reachable at this budget: %d of %d (%.0f%%)'
          % (reach, len(rows), 100.0*reach/len(rows)), flush=True)

    # the headline numbers
    Qall = max(r[2] for r in rows)
    Qmed = sorted(r[2] for r in rows)[len(rows)//2]
    print('\nlargest prime power anywhere in the 927 moduli : %d  -> k_max = %d'
          % (Qall, kmax(Qall, budget)), flush=True)
    print('median component largest prime power           : %d  -> k_max = %d'
          % (Qmed, kmax(Qmed, budget)), flush=True)
    print('\nTHE THRESHOLD, per component, as a table of Q against the largest affordable k:',
          flush=True)
    for k in range(2, 8):
        thr = int(budget ** (1.0/(k-1)))
        n = sum(1 for r in rows if r[2] <= thr)
        print('   k = %d  affordable while q^e <= %-10d  -> %d of %d components (%.0f%%)'
              % (k, thr, n, len(rows), 100.0*n/len(rows)), flush=True)

    json.dump({'budget': budget,
               'components': [{'n_cond': nc, 'n_wires': nw, 'largest_pp': Q, 'k_max': km}
                              for nc, nw, Q, km in rows]},
              open(os.path.join(V, 'v_cost.json'), 'w'), indent=1)
    print('\n-> v_cost.json', flush=True)
