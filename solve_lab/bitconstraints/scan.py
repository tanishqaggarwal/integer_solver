#!/usr/bin/env python3
"""Pass 2: inventory every atom / equation whose variable support lies inside
the 256 selector bits (or inside bits + exactly one extra variable)."""
import pickle, json, os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache.pkl')
CHAIN = os.path.join(HERE, '..', 'anneal', 'chain.json')


def load():
    with open(CACHE, 'rb') as f:
        D = pickle.load(f)
    chain = json.load(open(CHAIN))['chain_bit_vars']
    return D, chain, set(chain)


def support(poly):
    s = set()
    for m, c in poly:
        s.update(m)
    return s


def degree(poly):
    return max((len(m) for m, c in poly), default=0)


def fmt(poly):
    out = []
    for m, c in poly:
        if not m:
            out.append(f"{c:+d}")
        else:
            parts = []
            cnt = Counter(m)
            for v in sorted(cnt):
                parts.append(f"x_{v}" + (f"^{cnt[v]}" if cnt[v] > 1 else ""))
            body = "*".join(parts)
            if c == 1:
                out.append("+" + body)
            elif c == -1:
                out.append("-" + body)
            else:
                out.append(f"{c:+d}*{body}")
    s = "".join(out)
    return s[1:] if s.startswith('+') else s


def main():
    D, chain, BITS = load()
    atoms = D['atoms']
    eq_terms = D['eq_terms']
    eq_poly = D['eq_poly']
    eq_outer = D['eq_outer']
    idx = {v: i for i, v in enumerate(chain)}
    print(f"{len(atoms)} atoms, {len(eq_poly)} equations, {len(BITS)} bits")

    # --- where do bits appear at all ---
    atom_sup = [support(a) for a in atoms]
    atoms_touching = [i for i, s in enumerate(atom_sup) if s & BITS]
    print(f"atoms touching >=1 selector bit: {len(atoms_touching)}")

    bitonly_atoms = [i for i, s in enumerate(atom_sup) if s and s <= BITS]
    onestep_atoms = [i for i, s in enumerate(atom_sup)
                     if (s & BITS) and len(s - BITS) == 1]
    print(f"atoms with support SUBSET of bits : {len(bitonly_atoms)}")
    print(f"atoms with support = bits + 1 other: {len(onestep_atoms)}")

    # --- equations ---
    eq_sup = []
    for p in eq_poly:
        eq_sup.append(support(p))
    bitonly_eqs = [i for i, s in enumerate(eq_sup) if s and s <= BITS]
    onestep_eqs = [i for i, s in enumerate(eq_sup)
                   if (s & BITS) and len(s - BITS) == 1]
    print(f"EQUATIONS with support SUBSET of bits : {len(bitonly_eqs)}")
    print(f"EQUATIONS with support = bits + 1 other: {len(onestep_eqs)}")

    # atom -> equations containing it
    atom_eqs = defaultdict(list)
    for i, terms in enumerate(eq_terms):
        for c, aid in terms:
            atom_eqs[aid].append(i)

    # ---- report bit-only atoms grouped by shape ----
    shape = defaultdict(list)
    for aid in bitonly_atoms:
        a = atoms[aid]
        s = atom_sup[aid]
        d = degree(a)
        shape[(len(s), d, len(a))].append(aid)
    print("\nbit-only atom shapes  (nvars, degree, nterms) -> count")
    for k in sorted(shape):
        print(f"  {k} -> {len(shape[k])}")

    print("\n--- ALL bit-only atoms ---")
    for aid in sorted(bitonly_atoms,
                      key=lambda a: (len(atom_sup[a]), degree(atoms[a]))):
        eqs = atom_eqs[aid]
        print(f"atom#{aid} nv={len(atom_sup[aid])} deg={degree(atoms[aid])} "
              f"neq={len(eqs)} eqs={eqs[:8]}  {fmt(atoms[aid])[:300]}")

    print("\n--- bit-only EQUATIONS (hard constraints) ---")
    for i in bitonly_eqs:
        print(f"eq#{i} nv={len(eq_sup[i])} deg={degree(eq_poly[i])} "
              f"nterms={len(eq_poly[i])} outer={eq_outer[i]} "
              f"natoms={len(eq_terms[i])}")
        print("     ", fmt(eq_poly[i])[:600])

    json.dump({
        'bitonly_atoms': sorted(bitonly_atoms),
        'onestep_atoms': sorted(onestep_atoms),
        'bitonly_eqs': sorted(bitonly_eqs),
        'onestep_eqs': sorted(onestep_eqs),
    }, open(os.path.join(HERE, 'scan_index.json'), 'w'))
    print("\nwrote scan_index.json")


if __name__ == '__main__':
    main()
