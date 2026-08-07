"""Agent B: integer unit-propagation over the factored-atom model.

Goal: try to make EVERY atom zero (the natural 'circuit evaluates correctly' state),
seeded by the hard atoms (equations with |L|=1).
"""
import pickle, collections, sys, time

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model2.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']

fac_vars = []
for p in facs:
    s = set()
    for m in p: s.update(m)
    fac_vars.append(s)

# occurrence: var -> facs containing it
occ = collections.defaultdict(list)
for i, s in enumerate(fac_vars):
    for v in s: occ[v].append(i)

def evaluate(p, val):
    """substitute known vars; return (const, {var: linear_coef} or None if nonlinear-unknown)
       returns (k, lin, unknown_set, nonlinear_flag)"""
    k = 0
    lin = collections.defaultdict(int)
    unk = set()
    nonlin = False
    for m, c in p.items():
        cc = c
        rem = []
        for v in m:
            x = val.get(v)
            if x is None: rem.append(v)
            else: cc *= x
        if cc == 0:
            continue
        if not rem:
            k += cc
        elif len(rem) == 1:
            lin[rem[0]] += cc
            unk.add(rem[0])
        else:
            nonlin = True
            unk.update(rem)
    return k, lin, unk, nonlin

class Contradiction(Exception):
    pass

def propagate(target_facs, val, verbose=True):
    """target_facs: set of factor ids that must be 0. val: dict var->int (mutated)."""
    queue = collections.deque(target_facs)
    inq = set(target_facs)
    tf = set(target_facs)
    bad = set()
    t0 = time.time(); steps = 0
    while queue:
        f = queue.popleft(); inq.discard(f)
        steps += 1
        k, lin, unk, nonlin = evaluate(facs[f], val)
        if not unk:
            if k != 0:
                bad.add(f)
            continue
        if nonlin:
            continue
        if len(lin) == 1:
            v, c = next(iter(lin.items()))
            if c == 0:
                if k != 0: bad.add(f)
                continue
            if (-k) % c != 0:
                bad.add(f); continue
            x = (-k) // c
            if v in val:
                if val[v] != x: bad.add(f)
                continue
            val[v] = x
            for f2 in occ[v]:
                if f2 in tf and f2 not in inq:
                    queue.append(f2); inq.add(f2)
        if steps % 200000 == 0 and verbose:
            print("   steps=%d assigned=%d queue=%d %.1fs" % (steps, len(val), len(queue), time.time()-t0), flush=True)
    return bad

def main():
    # hard atoms
    hard_facs = set()
    for sc, L, k in eqs:
        if len(L) == 1 and sc != 0 and L[0][0] != 0:
            a = atoms[L[0][1]]
            if len(a) == 1:
                hard_facs.add(a[0])
    print("hard single-factor constraints:", len(hard_facs))
    val = {}
    bad = propagate(hard_facs, val)
    print("phase1: assigned=%d violated=%d" % (len(val), len(bad)))

    # phase 2: all single-factor atoms
    allf = set()
    for a in atoms:
        if len(a) == 1: allf.add(a[0])
    print("all single-factor constraints:", len(allf))
    val2 = dict(val)
    bad2 = propagate(allf, val2)
    print("phase2: assigned=%d violated=%d" % (len(val2), len(bad2)))
    pickle.dump({'val1': val, 'bad1': bad, 'val2': val2, 'bad2': bad2}, open(W+'prop.pkl','wb'), -1)
    # dump assignment json
    import json
    full = {('x_%d' % i): val2.get(i, 0) for i in range(38748)}
    json.dump(full, open(W+'out/prop_all0.json','w'))
    print("wrote out/prop_all0.json")

if __name__ == '__main__':
    main()
