"""S10 step 46: on the wire=1 branch, trade the two BIG checks for two 1-equation ones.

  x_11052 : d(a7930) = -1        d(a41512) = +27      a7930  is in 15 equations
  x_30163 : d(a29539) = -1       d(a40826) = +2       a29539 is in 12 equations

a40826 and a41512 each appear in exactly ONE equation, so zeroing a7930 and
a29539 at their expense is a huge net win in equation count.  Score by FAILING
EQUATIONS, not by number of nonzero atoms (that was the earlier greedy's mistake).
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

ROOTATOM = 37694
BLOCK = {ROOTATOM}


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            d = L.definer[u]
            if d in BLOCK:
                continue
            nv = T.solve_lin(d, u, v)
            if nv is not None:
                v[u] = nv
    return v


def rep(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'[{tag}] nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    for a in nz:
        print(f'     a{a:<6} neq={len(L.atom2eq.get(a,{})):<3} val={str(av[a])[:34]}')
    return av, nz, fail


v = L.load(os.path.join(HERE, 'wire1_solved2.json'))
av, nz, fail = rep(v, 'start')

for atom, knob in ((7930, 11052), (29539, 30163)):
    c, rest = T.lin_parts(atom, knob, v)
    if c and rest % c == 0:
        v[knob] = -rest // c
        fwd_block(v)
        av, nz, fail = rep(v, f'closed a{atom} via x_{knob}')
    else:
        print(f'  a{atom} not exactly solvable via x_{knob} (c={c})')

T.save(v, os.path.join(HERE, 'trade_out.json'))

print('\n=== what remains, and how expensive ===')
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
eqs = {}
for a in nz:
    eqs[a] = sorted(L.atom2eq.get(a, {}))
    print(f'  a{a}: {len(eqs[a])} equations {eqs[a][:14]}')
allq = set()
for a in nz:
    allq |= set(eqs[a])
print(f'  union: {len(allq)} equations')
fail = L.failing_eqs(av)
print(f'  actually failing: {len(fail)} -> score {L.NEQ-len(fail)}')
