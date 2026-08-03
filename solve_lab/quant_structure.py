#!/usr/bin/env python3
"""The twist is a DIOPHANTINE system between quantization levels:
  x_9770 = m*g,  x_3183 = m'*h        (22-side, m/m' small, from the codeword table)
  x_18274 = m2*g2, x_17728 = m2'*h2   (233-side)
Twist 1817 (with slack x_26977=x_20510*x_31302): 6033033*(m*g - m2*g2) + x_26977 = 0
Twist 44271 (RIGID): m'*h = m2'*h2
Compute all units, their gcds, and m-ranges to characterize/solve the Diophantine system."""
import json
from math import gcd
from functools import reduce
import numpy as np

def units(vals):
    vals = sorted(set(vals))
    G = reduce(gcd, [abs(v) for v in vals if v != 0]) if any(vals) else 0
    ms = sorted(set(v // G for v in vals)) if G else vals
    return G, ms

def main():
    cw = json.load(open('codewords.json'))
    x9770 = [int(x) for x in cw['x9770']]; x3183 = [int(x) for x in cw['x3183']]
    g, ms9770 = units(x9770)
    h, ms3183 = units(x3183)
    g2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    h2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    print(f"g  (x_9770 unit)  = {g}\n   bitlen {g.bit_length()}, m-values ({len(ms9770)}): {ms9770}")
    print(f"h  (x_3183 unit)  = {h}\n   bitlen {h.bit_length()}, m'-values ({len(ms3183)}): {ms3183}")
    print(f"g2 (x_18274 unit) = {g2}  bitlen {g2.bit_length()}")
    print(f"h2 (x_17728 unit) = {h2}  bitlen {h2.bit_length()}")
    print()
    print(f"gcd(g,  g2) = {gcd(g,g2)}  (bitlen {gcd(g,g2).bit_length()})")
    print(f"gcd(h,  h2) = {gcd(h,h2)}  (bitlen {gcd(h,h2).bit_length()})")
    print(f"gcd(g,  h)  = {gcd(g,h)}")
    print(f"gcd(g2, h2) = {gcd(g2,h2)}")
    print(f"gcd(g,  h2) = {gcd(g,h2)};  gcd(h, g2) = {gcd(h,g2)}")
    # Diophantine 44271: m'*h = m2'*h2.  gcd(h,h2)=d -> m' = (h2/d)*t, m2' = (h/d)*t
    d = gcd(h, h2)
    print(f"\n44271: m'*h = m2'*h2.  d=gcd(h,h2). m' must be multiple of h2/d = {h2//d} (bitlen {(h2//d).bit_length()})")
    print(f"  since m' in {ms3183[:5]}... (small), only m'=0 works unless h2/d is small.")
    # Diophantine 1817 with slack: 6033033*(m*g - m2*g2) = -x_26977 (a product value, ~free-ish)
    # This can bridge any m,m2 if x_26977 achievable. The RIGID one (44271) is the crux.
    print(f"\n=> If gcd-units force m'=0 (x_3183=x_17728=0), witness has those zero; 1817 slack bridges m*g vs m2*g2.")

if __name__ == '__main__':
    main()
