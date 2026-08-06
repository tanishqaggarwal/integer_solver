"""S11 step 70: the bipartite component of the seven failing equations.

Atoms and equations form a bipartite graph.  If the seven failing equations sit in
a SMALL component, the whole residual question is a self-contained subproblem that
can be settled exactly; if the component is the entire instance, it cannot.  Proper
BFS (the naive fixpoint loop is quadratic and never finished).
"""
import os, sys, json, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
FAIL = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
A, E = set(), set(FAIL)
q = collections.deque(('e', e) for e in FAIL)
while q:
    kind, x = q.popleft()
    if kind == 'e':
        for a in L.eq_atoms[x][2]:
            if a not in A:
                A.add(a)
                q.append(('a', a))
    else:
        for e in L.atom2eq[x]:
            if e not in E:
                E.add(e)
                q.append(('e', e))
print(f'component: {len(A)} atoms, {len(E)} equations '
      f'(instance has {L.NA} atoms, {L.NEQ} equations)')
print(f'  gate atoms {sum(1 for a in A if a in L.atom_out)}, '
      f'check atoms {sum(1 for a in A if a not in L.atom_out)}')
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
print(f'  nonzero atoms at the witness: {sum(1 for a in A if av[a])}')
json.dump({'atoms': sorted(A), 'eqs': sorted(E)},
          open(os.path.join(HERE, 'comp7.json'), 'w'))

# how far does the component grow if we stop after k hops?
A2, E2 = set(), set(FAIL)
front = list(FAIL)
for hop in range(8):
    na = set()
    for e in front:
        na |= set(L.eq_atoms[e][2])
    na -= A2
    A2 |= na
    ne = set()
    for a in na:
        ne |= set(L.atom2eq[a])
    ne -= E2
    E2 |= ne
    front = sorted(ne)
    print(f'  hop {hop+1}: {len(A2)} atoms, {len(E2)} equations')
    if not front:
        break
