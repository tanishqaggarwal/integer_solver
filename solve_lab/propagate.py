#!/usr/bin/env python3
"""Integer propagation engine over the atom set.

Loads all distinct atoms (polynomials that must equal 0).  Repeatedly:
  - substitute known variable values into each atom
  - if an atom reduces to a single unknown variable that is linear or a
    solvable quadratic, derive that variable's value (or finite domain)
  - propagate to fixpoint

Reports how many of the 38748 variables become determined.
"""
import json, time, sys
from collections import defaultdict

HERE = __file__.rsplit('/', 1)[0]
NVARS = 38748

def load_atoms():
    atoms = []
    with open(HERE + '/atoms/poly_atoms.jsonl') as f:
        for line in f:
            d = json.loads(line)
            poly = {tuple(m): c for m, c in d['poly']}
            atoms.append(poly)
    return atoms

def atom_vars(poly):
    s = set()
    for m in poly:
        s.update(m)
    return s

def substitute(poly, val):
    """Return reduced poly (dict monomial->coef) with known vars substituted."""
    out = defaultdict(int)
    for m, c in poly.items():
        newm = []
        coef = c
        ok = True
        for v in m:
            if val[v] is not None:
                coef *= val[v]
            else:
                newm.append(v)
        out[tuple(sorted(newm))] += coef
    return {m: c for m, c in out.items() if c != 0}

def solve_single(poly):
    """poly has exactly one unknown var. Return ('val', x) or ('dom', set) or
    ('contradiction', None) or ('skip', None)."""
    uvars = atom_vars(poly)
    assert len(uvars) == 1
    u = next(iter(uvars))
    # collect coefficients by power of u
    c0 = c1 = c2 = 0
    for m, c in poly.items():
        deg = len(m)  # all entries are u
        if deg == 0: c0 += c
        elif deg == 1: c1 += c
        elif deg == 2: c2 += c
        else:
            return ('skip', None)  # degree>2 in single var, skip
    if c2 == 0:
        if c1 == 0:
            return ('contradiction', None) if c0 != 0 else ('skip', None)
        if (-c0) % c1 != 0:
            return ('contradiction', None)
        return ('val', (u, (-c0) // c1))
    # quadratic c2 u^2 + c1 u + c0 = 0
    disc = c1 * c1 - 4 * c2 * c0
    if disc < 0:
        return ('contradiction', None)
    r = int(disc ** 0.5)
    while r * r > disc: r -= 1
    while (r + 1) * (r + 1) <= disc: r += 1
    if r * r != disc:
        return ('contradiction', None)  # no integer root
    roots = set()
    for s in (r, -r):
        num = -c1 + s
        if num % (2 * c2) == 0:
            roots.add(num // (2 * c2))
    if not roots:
        return ('contradiction', None)
    if len(roots) == 1:
        return ('val', (u, roots.pop()))
    return ('dom', (u, roots))

def main():
    t0 = time.time()
    atoms = load_atoms()
    print(f"loaded {len(atoms)} atoms in {time.time()-t0:.1f}s")
    var_atoms = defaultdict(list)
    for ai, poly in enumerate(atoms):
        for v in atom_vars(poly):
            var_atoms[v].append(ai)

    val = [None] * NVARS
    domain = {}   # var -> set of candidate ints (finite)
    contradictions = []

    # worklist: start with all atoms
    from collections import deque
    wl = deque(range(len(atoms)))
    inwl = [True] * len(atoms)

    def assign(v, x):
        if val[v] is not None:
            if val[v] != x:
                contradictions.append((v, val[v], x))
            return
        val[v] = x
        domain.pop(v, None)
        for ai in var_atoms[v]:
            if not inwl[ai]:
                inwl[ai] = True
                wl.append(ai)

    n_assigned_rounds = []
    steps = 0
    while wl:
        ai = wl.popleft()
        inwl[ai] = False
        poly = substitute(atoms[ai], val)
        uv = atom_vars(poly)
        steps += 1
        if len(uv) == 0:
            if poly.get((), 0) != 0:
                contradictions.append(('const', ai, poly))
            continue
        if len(uv) == 1:
            kind, data = solve_single(poly)
            if kind == 'val':
                u, x = data
                assign(u, x)
            elif kind == 'dom':
                u, roots = data
                if u in domain:
                    roots = domain[u] & roots
                domain[u] = roots
                if len(roots) == 1:
                    assign(u, next(iter(roots)))
            elif kind == 'contradiction':
                contradictions.append(('atom', ai, dict(poly)))
        # else: >=2 unknowns, wait

    n_assigned = sum(1 for x in val if x is not None)
    print(f"propagation fixpoint in {time.time()-t0:.1f}s, steps={steps}")
    print(f"assigned {n_assigned}/{NVARS} variables")
    print(f"vars with finite domain (unassigned): {len(domain)}")
    print(f"contradictions: {len(contradictions)}")
    for c in contradictions[:10]:
        print("   ", str(c)[:120])

    # save partial
    part = {f"x_{i}": val[i] for i in range(NVARS) if val[i] is not None}
    with open(HERE + '/partial_assignment.json', 'w') as g:
        json.dump(part, g)
    with open(HERE + '/domains.json', 'w') as g:
        json.dump({str(k): sorted(v) for k, v in domain.items()}, g)
    print(f"wrote partial_assignment.json ({len(part)} vars), domains.json")
    return val, domain, contradictions

if __name__ == '__main__':
    main()
