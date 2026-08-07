#!/usr/bin/env python3
"""agent V -- V1.  Reproduce L's |S|=17 single-wire run, then census the coupling structure.

Two things, both asked for and neither previously measured in general form:

(1) the SINGLE-WIRE end state at |S|=17, reported as NONZERO ATOMS OF 9,032 (never as a "stuck
    count", which drops entries whose residual is not 0 mod p);

(2) the connected components of the "shares a condition" graph, in BOTH readings, because they
    answer different questions:
       condition graph : nodes = undischarged c>1 conditions, edge iff they share a shift wire.
                         This is what L reported as [1,1].
       wire graph      : nodes = shift wires of those conditions, edge iff some condition is
                         influenced by both.  |component| = k is the exponent that prices the
                         joint solve: enumerate k-1 wires, root-find the last => q^(e(k-1)).

Usage: python3 v_diag.py <n | comma,list> [tag]
"""
import os, sys, json, time, collections, itertools
import v_base as B

V = '/home/user/integer_solver/solve_lab/agentV_work'
E, SL, p, SHIFT = B.E, B.SL, B.p, B.SHIFT


def single_wire_pass(vv, verbose=True):
    """L's closeS4 loop verbatim in semantics: per-wire group solve, global-nonzero guard."""
    gen = 0
    for outer in range(14):
        base = B.nzcount(vv)
        r = E.run(vv); gen += 1
        viol = B.violated(vv, r)
        if verbose:
            print('  outer %d: global nonzero %d, violated c-conditions %d'
                  % (outer, base, len(viol)), flush=True)
        if not viol:
            break
        wires = collections.defaultdict(list)
        for a in viol:
            for w in B.wires_of(a):
                wires[w].append(a)
        prog = 0
        for w, ats in sorted(wires.items(), key=lambda kv: -len(kv[1])):
            Vt = [a for a in ats if B.influences(vv, a, w)]
            if not Vt:
                continue
            t = B.GL['solve_group3'](vv, Vt, w, gen, base)
            if t:
                prog += 1
                base = B.nzcount(vv); gen += 1
        if verbose:
            print('     single-wire pass: %d accepted shift(s)' % prog, flush=True)
        if prog == 0:
            break
    return vv


def coupling(vv, atoms):
    """influence-checked wire sets, then components of both graphs."""
    A2W = {}
    for a in atoms:
        A2W[a] = sorted(w for w in B.wires_of(a) if B.influences(vv, a, w))
    W2A = collections.defaultdict(set)
    for a, ws in A2W.items():
        for w in ws:
            W2A[w].add(a)
    # union-find over conditions, joined by a shared wire
    par = {a: a for a in atoms}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x

    def uni(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            par[rx] = ry
    for w, ats in W2A.items():
        ats = sorted(ats)
        for b in ats[1:]:
            uni(ats[0], b)
    comps = collections.defaultdict(list)
    for a in atoms:
        comps[find(a)].append(a)
    out = []
    for root, ats in comps.items():
        ws = sorted(set(w for a in ats for w in A2W[a]))
        out.append({'conditions': ats, 'wires': ws})
    out.sort(key=lambda d: (-len(d['conditions']), -len(d['wires'])))
    return A2W, out


def report(vv, tag):
    # NOTE: solve_group3 restores vv[w] on a rejected trial but NOT the handle absorption that
    # nzcount's internal relift performed, so a final relift is REQUIRED before counting -- without
    # it the count reads 5 instead of 3 at |S|=17.  closeS4.close() does this; a caller that drives
    # the pass directly must too.
    B.relift(vv)
    r = E.run(vv)
    nz = [E.res[i] for i, x in enumerate(r) if x]
    viol = B.violated(vv, r)
    print('\n=== %s : NONZERO ATOMS = %d of 9032 ===' % (tag, len(nz)), flush=True)
    for a in nz:
        c = abs(SL[a])//p if SL.get(a) else 1
        print('    %-58s c=%d %s' % (a[:58], c,
              '(TARGET CONGRUENCE)' if any(t in a for t in B.TGTW) else ''), flush=True)
    print('  undischarged c>1 conditions: %d' % len(viol), flush=True)
    return nz, viol


if __name__ == '__main__':
    arg = sys.argv[1]
    S = [int(x) for x in arg.split(',')] if ',' in arg else B.onset(int(arg))
    tag = sys.argv[2] if len(sys.argv) > 2 else ('S%s' % arg)
    print('S (|S|=%d) = %s' % (len(S), S), flush=True)
    t0 = time.time()
    vv = B.greedy_init(S)
    print('greedy init done %.1f s -> global nonzero %d' % (time.time()-t0, B.nzcount(vv)), flush=True)
    single_wire_pass(vv)
    nz, viol = report(vv, '%s single-wire end state' % tag)

    # ---- coupling census on the SURVIVING conditions
    surv = [a for a in viol if not any(t in a for t in B.TGTW)]
    print('\n--- coupling census over the %d surviving non-target conditions ---' % len(surv),
          flush=True)
    A2W, comps = coupling(vv, surv)
    print('CONDITION-GRAPH COMPONENT SIZES: %s' % [len(d['conditions']) for d in comps], flush=True)
    print('WIRE-SET SIZE PER COMPONENT   : %s' % [len(d['wires']) for d in comps], flush=True)
    for i, d in enumerate(comps):
        print('  component %d: %d condition(s), %d wire(s)' % (i, len(d['conditions']), len(d['wires'])), flush=True)
        for a in d['conditions']:
            c = abs(SL[a])//p
            print('     c=%-12d %s' % (c, '*'.join('%d^%d' % qe for qe in sorted(B.factor(c).items()))), flush=True)
            print('        %s' % a[:96], flush=True)
        print('     wires: %s' % ['x%d' % w for w in d['wires']], flush=True)

    # ---- the same census over ALL 927 c>1 conditions (not just the survivors), which is the
    #      structure a general k-wire solver would have to range over
    print('\n--- coupling census over ALL c>1 conditions that any surviving wire touches ---',
          flush=True)
    allw = set(w for d in comps for w in d['wires'])
    nbr = [a for a in B.CGT2 if B.wires_of(a) & allw]
    A2W2, comps2 = coupling(vv, nbr)
    print('neighbourhood conditions: %d' % len(nbr), flush=True)
    print('CONDITION-GRAPH COMPONENT SIZES: %s' % [len(d['conditions']) for d in comps2], flush=True)
    print('WIRE-SET SIZE PER COMPONENT   : %s' % [len(d['wires']) for d in comps2], flush=True)

    json.dump({'S': S, 'nonzero_atoms': len(nz), 'nz': nz,
               'components': [{'conditions': d['conditions'], 'wires': d['wires']} for d in comps],
               'neighbourhood_components': [{'conditions': d['conditions'], 'wires': d['wires']}
                                            for d in comps2]},
              open(os.path.join(V, 'v_diag_%s.json' % tag), 'w'), indent=1)
    print('\nwall %.1f s -> v_diag_%s.json' % (time.time()-t0, tag), flush=True)
