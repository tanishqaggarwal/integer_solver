"""S10 step 79: the 1,800 variables in gate-DAG CYCLES -- a door never opened.

fwd.py reports the topological order covers only 29,675 of 31,475 defined
variables: 1,800 sit in cycles.  A cyclic block of gate equations is a SYSTEM,
not a chain, and a linear system with a nontrivial kernel has MANY solutions.
Every analysis so far picked one fixed point implicitly.  If a cycle admits a
kernel direction, that is free global movement no local method could see.

Find the strongly connected components of the gate graph, size them, and for each
non-trivial one test whether its local system is under-determined.
"""
import os, sys, collections, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
definer, atom_out = L.definer, L.atom_out

# gate graph: t -> depends on the other variables of its defining atom
dep = {}
for t, a in definer.items():
    dep[t] = [u for u in L.avars[a] if u != t and u in definer]

# Tarjan SCC (iterative)
index = {}; low = {}; onstack = {}; stack = []; comps = []; counter = [0]
for root in dep:
    if root in index:
        continue
    work = [(root, 0)]
    while work:
        node, pi = work[-1]
        if pi == 0:
            index[node] = low[node] = counter[0]; counter[0] += 1
            stack.append(node); onstack[node] = True
        recurse = False
        for i in range(pi, len(dep[node])):
            w = dep[node][i]
            if w not in index:
                work[-1] = (node, i + 1)
                work.append((w, 0)); recurse = True; break
            elif onstack.get(w):
                low[node] = min(low[node], index[w])
        if recurse:
            continue
        if low[node] == index[node]:
            comp = []
            while True:
                w = stack.pop(); onstack[w] = False; comp.append(w)
                if w == node: break
            comps.append(comp)
        work.pop()
        if work:
            parent = work[-1][0]
            low[parent] = min(low[parent], low[node])

big = [c for c in comps if len(c) > 1]
print(f'SCCs: {len(comps)} total; non-trivial (size > 1): {len(big)}')
sizes = collections.Counter(len(c) for c in big)
print(f'size histogram: {dict(sorted(sizes.items())[:15])}')
print(f'variables inside cycles: {sum(len(c) for c in big)}')

# For each non-trivial SCC, is the local gate system under-determined at v?
# Linearise: each t in the component has its gate atom; d(atom)/d(u) gives a matrix.
print('\n=== per-component local Jacobian rank (linearised at the witness) ===')
under = []
for comp in sorted(big, key=len)[:40]:
    cs = sorted(comp)
    idx = {u: i for i, u in enumerate(cs)}
    m = len(cs)
    J = [[0] * m for _ in range(m)]
    for i, t in enumerate(cs):
        a = definer[t]
        for u in L.avars[a]:
            if u in idx:
                # exact partial derivative at v
                s = 0
                for mo, c in L.polys[a].items():
                    k = mo.count(u)
                    if k == 0: continue
                    if k == 1:
                        term = c
                        for z in mo:
                            if z != u: term *= v[z]
                        s += term
                    else:
                        s += 2 * c * v[u]
                J[i][idx[u]] += s
    # rank over Q
    A = [[Fraction(x) for x in r] for r in J]
    r = 0
    for c in range(m):
        k = next((i for i in range(r, m) if A[i][c] != 0), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(m)]
        r += 1
        if r == m: break
    if r < m:
        under.append((m - r, m, cs))
    print(f'  size {m:>4}  rank {r:>4}  kernel {m-r:>3}'
          f'{"   <== UNDER-DETERMINED" if r < m else ""}')

print(f'\ncomponents with a nontrivial local kernel: {len(under)}')
for k, m, cs in under[:10]:
    print(f'   size {m}, kernel dim {k}, vars {cs[:10]}')
    # which atoms/equations does this component touch?
    ats = set()
    for t in cs:
        ats |= set(L.var_atoms[t])
    eqs = set()
    for a in ats:
        eqs |= set(L.atom2eq.get(a, ()))
    print(f'     touches {len(ats)} atoms, {len(eqs)} equations')
json.dump({'ncomps': len(big), 'sizes': dict(sizes),
           'under': [[k, m, cs[:20]] for k, m, cs in under]},
          open(os.path.join(HERE, 'cycles.json'), 'w'))
