"""Forward-evaluate the gate DAG from free inputs; report failing checks."""
import sys, collections, time
import dlib as L


def forward(v, block=()):
    """Evaluate every defined variable in topo order from the current free inputs.
    Returns number of gates that could not be solved exactly."""
    block = set(block)
    bad = 0
    for t in L.topo:
        a = L.definer[t]
        if a in block:
            continue
        nv = L.solve_for(a, t, v)
        if nv is None:
            bad += 1
            continue
        v[t] = nv
    return bad


def report(v, tag='', show=40):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    f = L.failing_eqs(av)
    print(f'[{tag}] score={L.NEQ-len(f)} failing={len(f)} nz_atoms={len(nz)}')
    ng = [a for a in nz if a in L.atom_out]
    nc = [a for a in nz if a not in L.atom_out]
    print(f'   nz gates={len(ng)} nz checks={len(nc)}')
    for a in nz[:show]:
        print(f'   a{a:<6} eqs={len(L.atom2eq.get(a,{})):<3} gate={a in L.atom_out} :: {L.atom_src[a][:160]}')
    return av, nz, f


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else None
    v = L.load(path) if path else [0] * L.NVARS
    # keep only free inputs
    for t in L.definer:
        v[t] = 0
    t0 = time.time()
    bad = forward(v)
    print(f'forward done in {time.time()-t0:.1f}s, unsolvable gates={bad}')
    report(v, tag=(path or 'zeros'))
