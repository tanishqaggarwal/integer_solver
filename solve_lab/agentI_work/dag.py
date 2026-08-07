#!/usr/bin/env python3
"""Build the atom/gate DAG: which atoms can define which variable, find free inputs."""
import pickle, os, collections, sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE, 'atoms.pkl'), 'rb'))
polys = pickle.load(open(os.path.join(HERE, 'polys.pkl'), 'rb'))
NV = 38748


def solvable_for(p):
    """Return set of vars v such that p is linear in v with coefficient +-1
    and v appears in no other monomial (so p=0 solves exactly for v over Z)."""
    occ = collections.defaultdict(list)
    for m, c in p.items():
        for v in set(m):
            occ[v].append((m, c))
    out = []
    for v, lst in occ.items():
        if len(lst) == 1:
            m, c = lst[0]
            if len(m) == 1 and abs(c) == 1:
                out.append(v)
    return out


sol = [solvable_for(p) for p in polys]
n_no = sum(1 for s in sol if not s)
print("atoms with no unit-solvable var:", n_no, "/", len(polys))

# candidate definers per variable
defs = collections.defaultdict(list)
for i, s in enumerate(sol):
    for v in s:
        defs[v].append(i)

vars_used = set()
for vs in D['atom_vars']:
    vars_used |= vs
print("vars appearing:", len(vars_used), "max id", max(vars_used))
print("vars with >=1 candidate definer:", len(defs))
print("vars with NO candidate definer (pure inputs):", len(vars_used - set(defs)))

cnt = collections.Counter(len(v) for v in defs.values())
print("definer-count histogram:", sorted(cnt.items())[:10])

pickle.dump({'sol': sol, 'defs': dict(defs), 'vars_used': vars_used},
            open(os.path.join(HERE, 'dag.pkl'), 'wb'))

# atoms with no solvable var -- what are they?
import re
sh = collections.Counter()
for i, s in enumerate(sol):
    if not s:
        t = re.sub(r'X\d+', 'V', D['atom_src'][i]); t = re.sub(r'\d+', 'N', t)
        sh[t] += 1
print("\nshapes of non-solvable atoms:")
for k, v in sh.most_common(20):
    print(f"  {v:6d}  {k}")
