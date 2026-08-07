#!/usr/bin/env python3
"""agent V -- V2.  The FULL component structure of the "shares a condition" graph.

L reported `[1,1]` at |S|=17.  That was the component structure of ONE residual pair at ONE point
in ONE search, not a statement about the instance.  This computes the whole thing, over all 927
c>1 conditions, in both readings, and prices the joint solve against it.

Two graphs, both needed and answering different questions:

  CONDITION graph : node = a c>1 condition, edge = they share a shift wire.
                    A component is a set of conditions that cannot be discharged independently.
  WIRE graph      : node = a shift wire that influences some c>1 condition, edge = some condition
                    is influenced by both.  |component| = k prices the joint solve at q^(e(k-1)).

Influence is PROBED (two direct recomputations per (atom,wire)), not assumed from the variable
sets -- a wire can appear in an atom's expression and have zero derivative.  Both the structural
and the probed graph are reported, because the difference is itself the measurement.

Usage: python3 v_comp.py <n | comma,list> [tag]
"""
import os, sys, json, time, collections, itertools
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, SHIFT = B.E, B.SL, B.p, B.SHIFT


def components(A2W):
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
    comps = collections.defaultdict(list)
    for a in A2W:
        comps[find(a)].append(a)
    out = []
    for ats in comps.values():
        ws = sorted(set(w for a in ats for w in A2W[a]))
        out.append((sorted(ats), ws))
    out.sort(key=lambda t: (-len(t[0]), -len(t[1])))
    return out


def summarise(name, comps):
    cs = collections.Counter(len(c) for c, w in comps)
    ws = collections.Counter(len(w) for c, w in comps)
    print('  %s: %d component(s)' % (name, len(comps)), flush=True)
    print('     condition-count histogram : %s' % dict(sorted(cs.items())), flush=True)
    print('     wire-count histogram      : %s' % dict(sorted(ws.items())), flush=True)
    print('     max conditions in one component = %d ; max wires = %d'
          % (max(len(c) for c, w in comps), max(len(w) for c, w in comps)), flush=True)
    return {'n_components': len(comps),
            'condition_hist': {str(k): v for k, v in sorted(cs.items())},
            'wire_hist': {str(k): v for k, v in sorted(ws.items())},
            'max_conditions': max(len(c) for c, w in comps),
            'max_wires': max(len(w) for c, w in comps)}


if __name__ == '__main__':
    arg = sys.argv[1]
    S = [int(x) for x in arg.split(',')] if ',' in arg else B.onset(int(arg))
    tag = sys.argv[2] if len(sys.argv) > 2 else ('S%s' % arg)
    print('S (|S|=%d) = %s' % (len(S), S), flush=True)
    t0 = time.time()
    vv = B.greedy_init(S)
    print('greedy fixpoint in %.1f s; global nonzero %d' % (time.time()-t0, B.nzcount(vv)), flush=True)
    r = E.run(vv)
    viol = B.violated(vv, r)
    print('violated c>1 conditions at the greedy fixpoint: %d of %d' % (len(viol), len(B.CGT2)),
          flush=True)

    res = {'S': S}

    # ---------- structural graph over ALL 927, no probing
    t1 = time.time()
    A2W_all = {a: sorted(B.wires_of(a)) for a in B.CGT2}
    print('\nSTRUCTURAL (all %d c>1 conditions, wires from the expression + value closure):'
          % len(B.CGT2), flush=True)
    ca = components(A2W_all)
    res['structural_all'] = summarise('all 927', ca)
    print('   built in %.1f s' % (time.time()-t1), flush=True)

    # ---------- probed graph over ALL 927
    t1 = time.time()
    A2W_pr = {}
    for a in B.CGT2:
        A2W_pr[a] = sorted(w for w in A2W_all[a] if B.influences(vv, a, w))
    A2W_pr = {a: ws for a, ws in A2W_pr.items() if ws}
    ndrop = len(B.CGT2) - len(A2W_pr)
    print('\nPROBED (influence tested by direct recomputation; %d condition(s) have NO influencing '
          'wire at all and drop out):' % ndrop, flush=True)
    cp = components(A2W_pr)
    res['probed_all'] = summarise('all with an influencing wire', cp)
    res['probed_no_wire'] = ndrop
    print('   built in %.1f s (%d probe pairs)' % (time.time()-t1, sum(len(v) for v in A2W_all.values())),
          flush=True)

    # ---------- the components the solver would actually have to range over
    print('\nCOMPONENTS CONTAINING A VIOLATED CONDITION (what a joint solve must handle):', flush=True)
    vs = set(viol)
    hot = [(c, w) for c, w in cp if vs & set(c)]
    if hot:
        res['hot'] = summarise('violated-carrying', hot)
        for c, w in hot:
            print('     %d condition(s) / %d wire(s); moduli %s'
                  % (len(c), len(w), sorted(abs(SL[a])//p for a in c))[:400], flush=True)
    else:
        print('     none (no violated condition has an influencing wire)', flush=True)

    # ---------- the cost table
    print('\nCOST: enumerate k-1 wires per prime power, root-find the last -> q^(e(k-1)) points.',
          flush=True)
    print('  largest prime power over all 927 moduli: ', end='', flush=True)
    big = 0
    for a in B.CGT2:
        for q, e in B.factor(abs(SL[a])//p).items():
            big = max(big, q**e)
    print('%d' % big, flush=True)
    res['largest_prime_power'] = big
    json.dump(res, open(os.path.join(V, 'v_comp_%s.json' % tag), 'w'), indent=1)
    print('\nwall %.1f s -> v_comp_%s.json' % (time.time()-t0, tag), flush=True)
