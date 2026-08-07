#!/usr/bin/env python3
"""Agent Z: complete classification of every ATOM that mentions a selector,
and of every EQUATION that mentions >=2 selectors."""
import os, sys, json, pickle, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly, pkey

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])

def classify(p, s):
    """p: polynomial dict; s: the (single) selector var in it. Return class tag."""
    mons = dict(p)
    lin = mons.get((s,), 0)
    sq = mons.get((s, s), 0)
    const = mons.get((), 0)
    others = {m: c for m, c in mons.items() if m not in ((s,), (s, s), ())}
    if not others:
        if sq != 0 and lin == -sq and const == 0:
            return 'BOOLEANITY'          # c*(s - s^2)
        if sq == 0 and lin != 0 and const == 0:
            return 'PIN0'                # c*s          (i.e. "s - 0")
        if sq == 0 and lin != 0 and const != 0:
            return 'PIN_CONST'           # c*s + d      (i.e. "s - 1")
        if sq != 0:
            return 'QUAD_OTHER'
        return 'CONST_ONLY'
    # atoms with other vars
    if len(others) == 1:
        (m, c), = others.items()
        if len(m) == 2 and s in m and lin == 0 and sq == 0 and const == 0:
            return 'PROD_SX'
        if len(m) == 2 and s in m:
            x = m[0] if m[1] == s else m[1]
            if mons.get((x,), 0) == -c and lin == 0 and sq == 0 and const == 0:
                return 'LOAD_OFF'        # x - s*x = (1-s)*x
    # s*(x-K)  ->  {(s,x):1, (s,):-K}
    if sq == 0 and const == 0 and len(others) == 1:
        (m, c), = others.items()
        if len(m) == 2 and s in m and lin != 0:
            return 'LOAD_ON'             # s*x - K*s = s*(x-K)
    return 'MIXED_OTHER'

def main():
    sys.setrecursionlimit(100000)
    lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]

    atom_seen = {}                    # pkey -> (poly, tag, selvars)
    tagcount = collections.Counter()
    sel_atom_tags = collections.defaultdict(collections.Counter)   # selector -> tags
    eq_sel_atoms = collections.defaultdict(list)   # eq -> [(coeff, pkey, tag, sel)]
    multi_sel_atoms = []

    for i, lhs in enumerate(lines):
        E = parse(lhs)
        L, _ = reduce_L(E)
        if not (varset(L) & sel):
            continue
        for c, a in atoms_of(L):
            p = poly(a)
            sv = set()
            for m in p:
                sv |= set(m) & sel
            if not sv:
                continue
            k = pkey(p)
            if len(sv) >= 2:
                multi_sel_atoms.append((i, c, k, sorted(sv)))
                tag = 'MULTISEL'
            else:
                s = next(iter(sv))
                tag = classify(p, s)
                sel_atom_tags[s][tag] += 1
            if k not in atom_seen:
                atom_seen[k] = (p, tag, sorted(sv))
                tagcount[tag] += 1
            eq_sel_atoms[i].append((c, k, tag, sorted(sv)))
        if i % 5000 == 0:
            print("  ...", i, flush=True)

    print()
    print("distinct atoms mentioning >=1 selector:", len(atom_seen))
    print("class histogram (distinct atoms):", dict(tagcount))
    print("atoms mentioning >=2 selectors:", len(multi_sel_atoms))
    print()
    # per selector
    percounts = collections.Counter()
    for s in sorted(sel):
        percounts[tuple(sorted(sel_atom_tags[s].items()))] += 1
    print("distinct per-selector class profiles:", len(percounts))
    for prof, n in percounts.most_common(20):
        print("   %4d selectors : %s" % (n, prof))

    # which selectors get a PIN atom
    pinned = [s for s in sel if sel_atom_tags[s].get('PIN0', 0) or sel_atom_tags[s].get('PIN_CONST', 0)]
    print()
    print("selectors appearing in a PIN0/PIN_CONST atom:", len(pinned))

    pickle.dump({'atom_seen': atom_seen, 'eq_sel_atoms': dict(eq_sel_atoms),
                 'sel_atom_tags': {k: dict(v) for k, v in sel_atom_tags.items()},
                 'multi_sel_atoms': multi_sel_atoms},
                open(os.path.join(HERE, 'zclass.pkl'), 'wb'))
    print("saved zclass.pkl")

if __name__ == '__main__':
    main()
