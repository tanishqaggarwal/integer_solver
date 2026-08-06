"""WR step 5: the exact linear model of a wire deformation.

Write w_u = p + d_u.  Every wire-identity atom (linear, all vars in the wire)
vanishes at d = 0, so its value IS its linear form in d.  Each equation carrying
such an atom therefore has value  form_e(d)  (or mult*form_e(d)^2 if squared) as
long as every non-identity atom in it stays zero.  Cost of a deformation d =
number of rows with form_e(d) != 0.
"""
import os, sys, collections, json, math
from fractions import Fraction
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
P = ad.P
ROOT, ROOTATOM = 26064, 37694

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
WIRE = W.wire_of(base)
widx = {u: i for i, u in enumerate(WIRE)}
N = len(WIRE)

IDENT = {}
NONLIN = []
for a in range(L.NA):
    vs = L.avars[a]
    if not vs or not all(u in widx for u in vs):
        continue
    form = collections.defaultdict(int)
    ok = True
    for m, c in L.polys[a].items():
        if len(m) == 1:
            form[widx[m[0]]] += c
        elif len(m) != 0:
            ok = False
            break
    if ok:
        IDENT[a] = {j: c for j, c in form.items() if c}
    else:
        NONLIN.append(a)

EQS = sorted(set().union(*[set(L.atom2eq.get(a, ())) for a in IDENT]))
rows = {}
for e in EQS:
    m, sq, co = L.eq_atoms[e]
    form = collections.defaultdict(int)
    for a, c in co.items():
        if a in IDENT:
            for j, cc in IDENT[a].items():
                form[j] += c * cc
    form = {j: c for j, c in form.items() if c}
    if form:
        rows[e] = form
RE = sorted(rows)


def cost_of(dvec, verbose=False):
    """dvec: list of length N. returns list of broken rows."""
    bad = []
    for e in RE:
        s = 0
        for j, c in rows[e].items():
            if dvec[j]:
                s += c * dvec[j]
        if s:
            bad.append(e)
    return bad


if __name__ == '__main__':
    print(f'wire {N}; identity atoms {len(IDENT)} (nonlinear wire-only atoms {len(NONLIN)})')
    print(f'equations touched {len(EQS)}, nontrivial rows {len(RE)}')
    print(f'a{ROOTATOM} in IDENT: {ROOTATOM in IDENT} -> {IDENT.get(ROOTATOM)}')
    print(f'a39417 in IDENT: {39417 in IDENT} -> {IDENT.get(39417)}')

    # rows whose ONLY support is the root coordinate
    rj = widx[ROOT]
    solo = [e for e in RE if set(rows[e]) == {rj}]
    print(f'rows supported on the ROOT coordinate alone: {solo}')

    # uniform deformation
    d = [1] * N
    bad = cost_of(d)
    print(f'\nUNIFORM d = (1,...,1): {len(bad)} rows break -> {bad}')

    # everything except the root
    d = [1] * N; d[rj] = 0
    bad2 = cost_of(d)
    print(f'ALL-BUT-ROOT           : {len(bad2)} rows break -> {bad2[:40]}')

    # root alone
    d = [0] * N; d[rj] = 1
    bad3 = cost_of(d)
    print(f'ROOT ALONE             : {len(bad3)} rows break -> {bad3}')

    # base failing equations of the deliverable -- overlap?
    F = W.F_WIRE
    b2 = list(base); F.fwd(b2)
    _, _, fail0, _ = F.report(b2, 'baseline', quiet=True)
    print(f'\nbaseline 7 failing eqs {sorted(fail0)}; intersect identity rows: '
          f'{sorted(set(fail0) & set(RE))}')

    json.dump({'wire': WIRE, 'rows': {str(e): {str(j): c for j, c in rows[e].items()}
                                      for e in RE}},
              open(os.path.join(HERE, 'wr_rows.json'), 'w'))
    print('saved wr_rows.json')
