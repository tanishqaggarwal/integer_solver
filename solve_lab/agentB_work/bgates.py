"""Agent B: gate/definition analysis on the factored-atom model (model2.pkl)."""
import pickle, collections, sys

M = pickle.load(open('/home/user/integer_solver/solve_lab/agentB_work/model2.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']

def fvars(p):
    s = set()
    for m in p:
        s.update(m)
    return s

def defines(p):
    """If poly p = c*(v - g(rest)) with v linear and not in any other monomial,
    return list of candidate output vars v (with coefficient)."""
    out = []
    for m, c in p.items():
        if len(m) != 1: continue
        v = m[0]
        # v must not appear in any other monomial
        ok = True
        for m2 in p:
            if m2 is m: continue
            if v in m2 and m2 != m:
                ok = False; break
        if ok:
            out.append((v, c))
    return out

def main():
    t = collections.Counter()
    # atom -> set of vars
    atom_vars = []
    for a in atoms:
        s = set()
        for f in a: s |= fvars(facs[f])
        atom_vars.append(s)

    # equations forcing a single atom to zero
    hard = collections.defaultdict(list)   # atom -> [eq idx]
    eq_atoms = []
    for e, (sc, L, k) in enumerate(eqs):
        eq_atoms.append([a for c, a in L])
        if len(L) == 1 and sc != 0 and L[0][0] != 0:
            hard[L[0][1]].append(e)
    print("equations with |L|=1:", sum(len(v) for v in hard.values()), " distinct atoms:", len(hard))

    # definition candidates from single-factor atoms
    defcand = collections.defaultdict(list)  # var -> [atom]
    single = 0
    for i, a in enumerate(atoms):
        if len(a) != 1: continue
        single += 1
        p = facs[a[0]]
        for v, c in defines(p):
            if abs(c) == 1:
                defcand[v].append(i)
    print("single-factor atoms:", single, " vars with >=1 def candidate:", len(defcand))
    h = collections.Counter(len(v) for v in defcand.values())
    print("def-candidate count hist:", sorted(h.items())[:12])

    allvars = set()
    for s in atom_vars: allvars |= s
    print("vars appearing:", len(allvars), "max", max(allvars))
    free = allvars - set(defcand)
    print("vars with NO def candidate (free inputs):", len(free))
    pickle.dump({'atom_vars': atom_vars, 'hard': dict(hard), 'defcand': dict(defcand),
                 'eq_atoms': eq_atoms},
                open('/home/user/integer_solver/solve_lab/agentB_work/gates.pkl','wb'), -1)
    print("wrote gates.pkl")

if __name__ == '__main__':
    main()
