"""Incremental forward evaluation on the CORRECTED engine (engine2)."""
import sys, os, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine2 as E2

NV = E2.NV; SEQ = E2.SEQ; definer = E2.definer; avars = E2.avars; acodes = E2.acodes
pos = {u: k for k, u in enumerate(SEQ)}

_UP = '/home/user/integer_solver/solve_lab/agentM_work/users2.pkl'
if os.path.exists(_UP):
    users = pickle.load(open(_UP, 'rb'))
else:
    users = collections.defaultdict(list)
    for w in SEQ:
        i, _ = definer[w]
        for u in avars[i]:
            if u != w:
                users[u].append(w)
    users = dict(users)
    pickle.dump(users, open(_UP, 'wb'))

atom_of = collections.defaultdict(list)
for i, vs in enumerate(avars):
    for u in vs:
        atom_of[u].append(i)


def downstream(changed):
    aff = set(); stack = list(changed)
    while stack:
        u = stack.pop()
        for w in users.get(u, ()):
            if w not in aff:
                aff.add(w); stack.append(w)
    return aff


def apply_delta(v0, changes):
    v = list(v0)
    for k, val in changes.items():
        v[k] = val
    aff = downstream(changes.keys())
    ns = {'v': v, '__builtins__': {}}
    for u in sorted(aff, key=lambda u: pos[u]):
        i, kind = definer[u]
        E2._solvevar(v, ns, u, i, kind[0])
    return v, aff


def atoms_touching(aff):
    s = set()
    for u in aff:
        s.update(atom_of[u])
    return s


def resid_delta(v0, base_bad, changes):
    """Return (new bad-atom dict, new vector)."""
    v, aff = apply_delta(v0, changes)
    touched = atoms_touching(set(aff) | set(changes))
    ns = {'v': v, '__builtins__': {}}
    bad = dict(base_bad)
    for i in touched:
        r = eval(acodes[i], ns)
        if r:
            bad[i] = r
        else:
            bad.pop(i, None)
    return bad, v
