#!/usr/bin/env python3
"""Pass 15: how much booleanity of the selectors is a RIGOROUS consequence of
whole equations (no "every atom vanishes" assumption)?

For integer x,  t(x) = x - x^2 = x(1-x) <= 0,  with equality iff x in {0,1}.
So an equation  sum_i c_i * t(x_i) = 0  with all c_i of the same sign forces
t(x_i) = 0 for every i, i.e. every x_i is boolean.
Once some variables are known boolean their terms vanish, which can make other
equations sign-uniform -- so this propagates.
"""
import os, json
from collections import defaultdict
from scan import load, support, degree

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
eq = D['eq_poly']
idx = json.load(open(os.path.join(HERE, 'scan_index.json')))


def boolvec(poly):
    d = defaultdict(lambda: [0, 0])
    for m, c in poly:
        if len(m) == 1:
            d[m[0]][0] += c
        elif len(m) == 2 and m[0] == m[1]:
            d[m[0]][1] += c
        else:
            return None
    out = {}
    for v, (a, b) in d.items():
        if a != -b:
            return None
        if a:
            out[v] = a
    return out


rows = []
for e in range(len(eq)):
    bv = boolvec(eq[e])
    if bv:
        rows.append((e, bv))
print(f"equations that are pure Z-combinations of booleanity terms t(x)=x-x^2: "
      f"{len(rows)}")
cov = set()
for e, bv in rows:
    cov |= set(bv)
print(f"variables covered: {len(cov)} ({len(cov & BITS)} of them selectors)")

known = set()
prooflog = {}
changed = True
while changed:
    changed = False
    for e, bv in rows:
        rem = {v: c for v, c in bv.items() if v not in known}
        if not rem:
            continue
        signs = {1 if c > 0 else -1 for c in rem.values()}
        if len(signs) == 1:
            for v in rem:
                known.add(v)
                prooflog[v] = e
            changed = True
print(f"\nvariables PROVEN boolean by sign-uniform integer equations: {len(known)}")
print(f"   of which selector bits: {len(known & BITS)} / 256")
by_eq = defaultdict(list)
for v, e in prooflog.items():
    by_eq[e].append(v)
print(f"   certificates: {len(by_eq)} equations")
for e in sorted(by_eq)[:12]:
    print(f"     eq#{e} -> proves booleanity of "
          f"{sorted('x_%d' % v for v in by_eq[e])}")

json.dump({'n_bool_combination_equations': len(rows),
           'covered_vars': sorted(cov),
           'proven_boolean': sorted(known),
           'proven_boolean_selectors': sorted(known & BITS),
           'certificate_equation_per_var': {str(k): v for k, v in prooflog.items()}},
          open(os.path.join(HERE, 'booleanity_proof.json'), 'w'))
print("wrote booleanity_proof.json")
