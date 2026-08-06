"""S10 step 49: deviate only the wire members that actually matter.

The whole-wire-at-1 branch costs 13 equations because it breaks the root pin
a37694.  But I only need the handles used by a7930 and a29539 to be unquantised,
i.e. only these two wire members need to leave p:

    a7929  = x_7927  - x_15616 * x_11052     (handle of a7930)
    a29538 = x_29967 - x_11360 * x_30163     (handle of a29539)

Deviating a single member breaks only the copy atom that defines it, leaving the
root pin intact.  Measure the cost of each option.
"""
import os, sys, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIRE = set(u for u in range(L.NVARS) if base[u] == P)
print(f'wire size {len(WIRE)}')

for u in (15616, 11360):
    d = L.definer.get(u)
    print(f'\nx_{u}: definer a{d} = {L.atom_src[d]}   neq={len(L.atom2eq.get(d,{}))}')
    print(f'   occurs in {len(L.var_atoms[u])} atoms, {len(L.var_eqs[u])} equations')
    # who else copies FROM u?
    users = [a for a in L.var_atoms[u]
             if a in L.atom_out and L.atom_out[a][1] != u
             and all(w in WIRE for w in L.avars[a])]
    print(f'   wire members copied FROM x_{u}: {[L.atom_out[a][1] for a in users]}')


def fwd_block(v, block, rounds=4):
    for _ in range(rounds):
        for x in ad.ORDER:
            a = L.definer[x]
            if a in block:
                continue
            nv = T.solve_lin(a, x, v)
            if nv is not None:
                v[x] = nv
    return v


def rep(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'[{tag}] nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}', flush=True)
    return av, nz, fail


# --- option A: deviate only x_15616 and x_11360 (and their copy-descendants) --
def descendants(u):
    out, frontier = {u}, [u]
    while frontier:
        nxt = []
        for w in frontier:
            for a in L.var_atoms[w]:
                if a not in L.atom_out:
                    continue
                t = L.atom_out[a][1]
                if t != w and t in WIRE and t not in out and all(z in WIRE for z in L.avars[a]):
                    out.add(t); nxt.append(t)
        frontier = nxt
    return out


for tag, seeds in (('x_15616 only', {15616: 1}),
                   ('x_11360 only', {11360: 1}),
                   ('both', {15616: 1, 11360: 1})):
    v = list(base)
    block = set()
    touched = set()
    for u, val in seeds.items():
        dset = descendants(u)
        touched |= dset
        for t in dset:
            v[t] = val
            block.add(L.definer[t])
    fwd_block(v, block)
    print(f'\n--- {tag}: deviated {len(touched)} wire members, blocked {len(block)} copy atoms')
    rep(v, tag)
    T.save(v, os.path.join(HERE, f'partial_{tag.split()[0].replace("x_","")}.json'))
