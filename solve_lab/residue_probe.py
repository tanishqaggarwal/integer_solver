#!/usr/bin/env python3
"""Cheap structural probe of the target values and residue pools:
GCDs, shared factors, small-modulus congruences. Looking for a hidden modulus or
a lattice relation connecting the 22-side and 233-side."""
import json
from math import gcd
from functools import reduce

T = {
 '9770':  119182891324903069288022589460020572593207162963685444009526935473255725746139626528632451,
 '3183':  62388561396646277577754745632050284891732896280504135411098502190267395125606332503618812,
 '18274': 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002,
 '17728': 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626,
}

def main():
    rp = json.load(open('residue_pool.json'))
    coarse = [int(x) for x in rp['pool22']]
    fine = [int(x) for x in rp['pool233']]
    print(f"coarse residues: {len(coarse)}, fine: {len(fine)}")

    D = T['9770'] - T['18274']     # twist gap 1 (must be absorbed)
    D2 = T['3183'] - T['17728']    # twist gap 2
    print(f"gap D1 = x_9770-x_18274 = {D}")
    print(f"gap D2 = x_3183-x_17728 = {D2}")
    print(f"bitlen D1={D.bit_length()} D2={D2.bit_length()}")

    gall = reduce(gcd, coarse + fine)
    print(f"gcd(all residues) = {gall}")
    print(f"gcd(D1, all residues) = {reduce(gcd, coarse+fine, abs(D))}")
    print(f"gcd(D1, D2) = {gcd(abs(D), abs(D2))}")
    # pairwise gcd stats of residues (are they coprime or share factors?)
    import itertools
    sh = 0; tot = 0
    for a, b in itertools.islice(itertools.combinations(coarse+fine, 2), 3000):
        g = gcd(a, b); tot += 1
        if g > 1: sh += 1
    print(f"residue pairs (sampled {tot}) sharing a common factor >1: {sh}")

    # small-prime congruence signature of the targets and D
    smallp = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]
    for name, v in list(T.items()) + [('D1', D), ('D2', D2)]:
        sig = [v % p for p in smallp]
        print(f"  {name} mod small primes: {sig}")

    # is D1 a small multiple of any single residue? (D1 / r near integer)
    hits = [(r, D % r) for r in coarse+fine if D % r == 0]
    print(f"residues dividing D1 exactly: {len(hits)}")
    # gcd(D1, product-ish): check gcd(D1, each residue) large?
    big_g = sorted(((gcd(abs(D), r), r) for r in coarse+fine), reverse=True)[:3]
    print(f"largest gcd(D1, residue): {[(g, str(r)[:20]+'...') for g,r in big_g]}")

if __name__ == '__main__':
    main()
