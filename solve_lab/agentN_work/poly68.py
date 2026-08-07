"""Redo the knob SELECTION symbolically too — the last place a finite difference was still used.

`widen.wide_knobs` starts from the exact syntactic support (68 free inputs) but then keeps only
those whose step-1 bump moves an atom of the region.  That is a finite-difference filter: a free
input whose first-order effect cancels against its second-order effect at t = 1 is dropped even
though it genuinely moves the region.  Here the filter is replaced by the exact polynomial support:
a free input is a knob iff it appears in the exact polynomial of some region row.

Then the whole reduction is redone on that complete knob set.
"""
import os, sys, json, time, pickle
from collections import defaultdict
from flint import fmpz_mat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN, zsolve
from optN import make, build, WIT, fr, FREE, FR0, atom_eqs, _bits, inner
import frameB as FB
from polyexact import P
from polyfull import exact_polys, evalP, parts
from kerquad import int_kernel_columns

eq_terms = ev.eq_terms


def main():
    st = make(WIT)
    b0 = build(st)
    Rl = b0['R']
    Rset = set(Rl)

    # exact syntactic candidate set (same starting point widen.py uses)
    atoms_R = set()
    for e in Rl:
        for c, a in eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    cands = sorted(y for y in cands if y in FREE)
    print('exact syntactic candidates: %d' % len(cands), flush=True)

    # every equation any candidate can touch, by exact syntactic support
    touched = set()
    for Y in cands:
        for a in fr.chk.get(Y, []):
            touched |= atom_eqs[a]
    outside = sorted(touched - Rset)
    print('equations any candidate can touch: %d region + %d outside'
          % (len(Rl), len(outside)), flush=True)

    rows = list(Rl) + list(outside)
    t0 = time.time()
    polys = exact_polys(st, rows, cands)
    print('exact polynomials over all %d candidates built in %.1fs' % (len(cands), time.time() - t0),
          flush=True)

    k = len(cands)
    PARTS = {}
    degdist = defaultdict(int)
    for e in rows:
        c0 = 0
        lin = [0] * k
        hi = {}
        for mono, v in polys[e].c.items():
            dg = sum(mono)
            if dg == 0:
                c0 = v
            elif dg == 1:
                lin[mono.index(1)] = v
            else:
                hi[mono] = v
        PARTS[e] = (c0, lin, hi)
        degdist[polys[e].deg()] += 1
    print('\nDEGREE DISTRIBUTION over all %d rows: %s'
          % (len(rows), dict(sorted(degdist.items()))), flush=True)
    print('max total degree of the exact system: %d' % max(degdist), flush=True)
    rdeg = {e: polys[e].deg() for e in Rl}
    print('region row degrees: %s' % rdeg, flush=True)
    print('max terms in a row: %d ; max coefficient bits: %d'
          % (max(polys[e].nterms() for e in rows),
             max(polys[e].maxbits() for e in rows)), flush=True)

    # which candidates genuinely appear in a region row's polynomial?
    used = set()
    for e in Rl:
        c0, lin, quad = PARTS[e]
        for j, c in enumerate(lin):
            if c:
                used.add(j)
        for mono in quad:
            for j, x in enumerate(mono):
                if x:
                    used.add(j)
    print('\ncandidates that genuinely move the region (exact polynomial support): %d of %d'
          % (len(used), k), flush=True)
    widen_set = set(json.load(open(os.path.join(HERE, 'runs', 'polyfull.json')))['knobs'])
    extra = sorted(cands[j] for j in used if cands[j] not in widen_set)
    lost = sorted(y for y in widen_set if cands.index(y) not in used)
    print('   knobs the step-1 filter MISSED: %d -> %s' % (len(extra), extra), flush=True)
    print('   knobs the step-1 filter kept that do NOT move the region: %d -> %s'
          % (len(lost), lost), flush=True)

    # ---- redo the reduction on the complete set ------------------------------------------
    print('\n=== reduction on the complete knob set (%d unknowns) ===' % k, flush=True)
    outside_live = [e for e in outside if any(PARTS[e][1]) or PARTS[e][2] or PARTS[e][0]]
    print('outside rows the candidates actually move: %d' % len(outside_live), flush=True)
    assert all(PARTS[e][0] == 0 for e in outside_live), 'outside row nonzero at the witness'

    LIN = [e for e in outside_live if not PARTS[e][2]]
    QUA = [e for e in outside_live if PARTS[e][2]]
    print('   %d purely linear, %d quadratic' % (len(LIN), len(QUA)), flush=True)

    K0 = int_kernel_columns([PARTS[e][1] for e in LIN], k)
    Kr = fmpz_mat([[int(x) for x in v] for v in K0]).lll().tolist()
    K0 = [[int(x) for x in r] for r in Kr if any(r)]
    d = len(K0)
    print('ker_Z of the linear collateral rows: dim %d' % d, flush=True)

    def restrict(e, K):
        """substitute t = K s exactly, at any degree"""
        dd = len(K)
        P.NK = dd
        T = [P({tuple(1 if a == b else 0 for b in range(dd)): K[a][j]
                for a in range(dd) if K[a][j]}) for j in range(k)]
        c0, lin, hi = PARTS[e]
        acc = P.const(c0)
        for j, c in enumerate(lin):
            if c:
                acc = acc + T[j] * c
        for mono, c in hi.items():
            term = P.const(c)
            for j, ex in enumerate(mono):
                for _ in range(ex):
                    term = term * T[j]
            acc = acc + term
        L = [0] * dd
        Q = {}
        cc = 0
        for mono, v in acc.c.items():
            g = sum(mono)
            if g == 0:
                cc = v
            elif g == 1:
                L[mono.index(1)] = v
            else:
                idx = [a for a, x in enumerate(mono) if x]
                key = (idx[0], idx[0]) if len(idx) == 1 and mono[idx[0]] == 2 else tuple(
                    a for a, x in enumerate(mono) for _ in range(x))
                Q[key] = v
        return cc, L, Q

    forced = set()
    branch = []
    for e in QUA:
        c0, L, Q = restrict(e, K0)
        nzl = [a for a, c in enumerate(L) if c]
        if not nzl and not Q:
            continue
        if c0 == 0 and len(nzl) == 1 and not Q:
            forced.add(nzl[0])
        elif c0 == 0 and not nzl and len(Q) == 1 and list(Q)[0][0] == list(Q)[0][1]:
            forced.add(list(Q)[0][0])
        else:
            branch.append((e, c0, L, Q))
    print('quadratic collateral rows that force a coordinate to 0: %d (coords %s)'
          % (len(forced), sorted(forced)), flush=True)
    print('quadratic collateral rows that need a genuine branch: %d' % len(branch), flush=True)
    for e, c0, L, Q in branch:
        print('   eq %-6d const=%s lin=%d quad=%d %s'
              % (e, c0 != 0, sum(1 for x in L if x), len(Q), sorted(Q)[:6]), flush=True)

    free = [a for a in range(d) if a not in forced]
    M, b, nq = [], [], 0
    for e in Rl:
        c0, L, Q = restrict(e, K0)
        if Q:
            nq += 1
        M.append([L[a] for a in free])
        b.append(c0)
    print('region rows genuinely quadratic on the variety: %d of 12' % nq, flush=True)
    if nq == 0 and not branch:
        opt, rws, exh, nd = zsolve.max_zero_rows(M, b, len(free), len(M), node_cap=4000000)
        fail = len(M) - opt
        print('\nOPT = %d of %d on a rank-%d lattice, exhaustive=%s, nodes=%d'
              % (opt, len(M), len(free), exh, nd), flush=True)
        print('   rows: %s' % [Rl[i] for i in rws], flush=True)
        print('   => failing = %d, score = %d' % (fail, 39033 - fail), flush=True)
        json.dump(dict(cands=len(cands), used=len(used), extra=extra, lost=lost,
                       kernel=d, forced=sorted(forced), rank=len(free), opt=opt,
                       exhaustive=bool(exh), failing=fail, score=39033 - fail),
                  open(os.path.join(HERE, 'runs', 'poly68.json'), 'w'), indent=1)
    else:
        json.dump(dict(cands=len(cands), used=len(used), extra=extra, lost=lost,
                       kernel=d, forced=sorted(forced), branches=len(branch), region_quad=nq),
                  open(os.path.join(HERE, 'runs', 'poly68.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
