#!/usr/bin/env python3
"""Extract the huge-atom (residue-load) network directly from the atom set,
independent of any forward-eval orientation. huge-atom form:
    bit*(x_B - HUGE) - s*x_C = 0   (bit in {0,1}) i.e.  bit*x_B - bit*HUGE - s*x_C
We want the TRUE algebraic assembly of x_18274 and x_17728 from the 233 bits.
"""
import json
from collections import defaultdict
from confluent_eval5 import boolean_vars
from propagate import load_atoms, atom_vars

def main():
    A = load_atoms()
    bset = boolean_vars(A)
    control = set(json.load(open('control_bits.json')))
    BITS22 = {1782,1858,2795,2800,3483,5443,10652,19520,21188,21588,23634,26947,27512,29682,30104,30596,30658,30792,33251,37748,37885,38116}

    # classify huge atoms: exactly monomials with a big (>1e40) constant
    huge = []
    for a, poly in enumerate(A):
        bigterms = [(m, c) for m, c in poly.items() if abs(c) > 10**40]
        if not bigterms: continue
        huge.append(a)
    print(f"atoms with a huge (>1e40) coefficient: {len(huge)}")

    # for each, describe: degree-2 terms (bit*x_B), the huge const term, linear terms
    forms = []
    control_bits_used = set()
    for a in huge:
        poly = A[a]
        prods = [(m, c) for m, c in poly.items() if len(m) == 2]
        lins = [(m, c) for m, c in poly.items() if len(m) == 1]
        consts = [(m, c) for m, c in poly.items() if len(m) == 0]
        bigconst = [(m, c) for m, c in lins + prods if abs(c) > 10**40]
        forms.append((a, len(prods), len(lins), len(consts), len(poly)))
    # histogram of shapes
    from collections import Counter
    shape = Counter((p, l, c) for a, p, l, c, n in forms)
    print("shape (nprod,nlin,nconst) histogram (top 15):")
    for sh, cnt in shape.most_common(15):
        print(f"   nprod={sh[0]} nlin={sh[1]} nconst={sh[2]}: {cnt}")

    # Which huge atoms mention a control bit (as a var)?
    hb = 0
    bit_to_huge = defaultdict(list)
    for a in huge:
        va = atom_vars(A[a])
        cb = va & control
        if cb:
            hb += 1
            for b in cb: bit_to_huge[b].append(a)
    print(f"huge atoms mentioning a control bit: {hb}")
    print(f"control bits appearing in huge atoms: {len(bit_to_huge)} / {len(control)}")

    # sample a few huge atoms fully
    print("\nsample huge atoms:")
    for a in huge[:4]:
        poly = A[a]
        print(f" atom {a}: {len(poly)} terms")
        for m, c in sorted(poly.items(), key=lambda kv: (-len(kv[0]), kv[0])):
            cb = [x for x in m if x in control]
            tag = f"  <bit {cb}>" if cb else ""
            cs = str(c) if abs(c) < 10**12 else f"{'HUGE' if abs(c)>10**40 else 'big'}({len(str(abs(c)))}d)"
            vv = '*'.join(f'x_{x}' for x in m) if m else '1'
            print(f"     {cs} * {vv}{tag}")

    # Do x_18274 / x_17728 appear directly in any huge atom?
    for w in (18274, 17728, 9770, 3183):
        ha = [a for a in huge if w in atom_vars(A[a])]
        print(f"x_{w} appears in {len(ha)} huge atoms: {ha[:6]}")

    json.dump({'huge_atoms': huge,
               'bit_to_huge': {str(b): v for b, v in bit_to_huge.items()}},
              open('huge_network.json', 'w'))
    print("wrote huge_network.json")

if __name__ == '__main__':
    main()
