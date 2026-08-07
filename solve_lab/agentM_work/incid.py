"""How thin is the incident set?  Characterised from the EQUATION side only.

I do not enumerate sites -- that stays L's job.  Instead: a handle can only help if moving
it changes some atom that appears in a baseline-failing equation.  Moving a freed handle u
changes EVERY atom that mentions u (occ[u]), so

    u is incident  <=>  occ[u]  intersects  {atoms appearing in the baseline failures}

That is a necessary condition, computable for every variable at once, with no site
enumeration.  Its size is the answer to "is incidence forced onto the deliverable's site or
is it a property many sites could have".
"""
import sys, os, re, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import price as PR

BS = json.load(open('baseline_sets.json'))
FA, FB = BS['A'], BS['B']
INTER = BS['intersection']
USE = INTER if not BS['identical'] else FA
print(f'baselines identical: {BS["identical"]}; using '
      f'{"INTERSECTION" if not BS["identical"] else "the common set"} of {len(USE)} equations')

# atoms appearing in those equations (with nonzero coefficient)
A = set()
percount = collections.Counter()
for e in USE:
    for c, a in H.eqt[e][2]:
        if a >= 0 and c:
            A.add(a)
            percount[a] += 1
print(f'\natoms appearing in the baseline-failing equations: {len(A)}')

# which variables can move them?
inc = collections.defaultdict(set)
for a in A:
    for u in H.avars[a]:
        inc[u].add(a)
print(f'variables that touch at least one such atom: {len(inc)}')

# restrict to plausible HANDLES: product-defined variables (h = a*b), P's population shape
prod = set()
bare = set()
for u in H.SEQ:
    a = H.definer[u][0]
    t = H.atoms[a]
    if re.fullmatch(r'x_%d - x_\d+ \* x_\d+' % u, t):
        prod.add(u)
    elif t == 'x_%d' % u:
        bare.add(u)
print(f'product-defined variables in my frame: {len(prod)}   bare-defined: {len(bare)}')

cand = (prod | bare)
incident_handles = sorted(u for u in cand if u in inc)
print(f'\n*** INCIDENT HANDLE POOL: {len(incident_handles)} of {len(cand)} '
      f'({100.0*len(incident_handles)/max(1,len(cand)):.2f}%) ***')

D4 = [642, 28730, 29854, 31864]
print(f'the deliverable\'s four in the pool: { {u: (u in inc) for u in D4} }')
print(f'the deliverable\'s four are product/bare-defined: { {u: (u in cand) for u in D4} }')

rank = sorted(incident_handles, key=lambda u: -len(inc[u]))
print(f'\nincident handles ranked by how many failing-equation atoms they touch:')
for u in rank[:40]:
    mark = '  <== DELIVERABLE' if u in D4 else ''
    neq = len(set(e for e in USE
                  for c, a in H.eqt[e][2] if a >= 0 and c and a in inc[u]))
    print(f'  x_{u:<6d} touches {len(inc[u]):2d} atoms, {neq:2d} of the {len(USE)} equations{mark}')

json.dump({'n_equations': len(USE), 'equations': USE,
           'n_atoms': len(A), 'n_candidate_handles': len(cand),
           'incident_handles': incident_handles,
           'deliverable_four': D4},
          open('incident_pool.json', 'w'), indent=1)
print(f'\nwrote incident_pool.json ({len(incident_handles)} incident handles)')
