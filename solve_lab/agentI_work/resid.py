import sys, os
from model import Model, load_assign
M = Model()
v = load_assign(sys.argv[1])
fails, av, cv = M.eq_fail(v)
nz = [i for i, x in enumerate(av) if x != 0]
print("nonzero atoms:")
for i in nz:
    print(f"  a{i}: {M.src[i]}  = {av[i]}")
    print(f"        in {len(M.atom_eqs[i])} eqs")
print("\nfailing equations, their atom decomposition:")
for e in fails:
    print(f"  eq {e}: outer={M.eq_outer[e]}  core={cv[e]}")
    for c, a in M.eq_terms[e]:
        if av[a]:
            print(f"      {c:+d} * a{a}[{M.src[a]}] = {av[a]}")
print("\nequations touched by nonzero atoms (all):")
touched = set()
for i in nz:
    for e, c in M.atom_eqs[i]:
        touched.add(e)
print(sorted(touched), len(touched))
for e in sorted(touched):
    parts = [(c, a, av[a]) for c, a in M.eq_terms[e] if av[a]]
    print(f"  eq{e} core={cv[e]} :", [(c, a) for c, a, _ in parts])
