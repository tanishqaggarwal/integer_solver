#!/usr/bin/env python3
"""CREATIVE data-mining of the construction constants. The setter built atoms that
vanish at the witness; the HUGE residue-load constants and the base quanta g,g2,h,h2
encode it. Hunt for exploitable structure:
 - factor g,g2,h,h2 (trial division + Pollard rho) -> small/shared factors?
 - relations among g,g2,h,h2 (h2/h, g/g2, integer combos)?
 - the HUGE load constants: gcds with g2/h2, ratios, and whether the twist quantum
   appears as a simple combination.
"""
import json, time, math
from collections import Counter
from propagate import load_atoms, atom_vars

def rho(n):
    if n % 2 == 0: return 2
    x = 2; y = 2; d = 1; c = 1
    f = lambda v: (v*v + c) % n
    while d == 1:
        x = f(x); y = f(f(y)); d = math.gcd(abs(x-y), n)
    return d if d != n else None

def factor(n, budget=200000):
    n = abs(n); fs = Counter()
    for p in range(2, budget):
        while n % p == 0: fs[p] += 1; n //= p
        if p*p > n: break
    if n > 1:
        # a few rho attempts for remaining
        stack = [n]; tries = 0
        while stack and tries < 40:
            m = stack.pop(); tries += 1
            if m == 1: continue
            # simple primality-ish: if small enough already handled
            d = rho(m) if m > 1 else None
            if d and d != m:
                stack.append(d); stack.append(m//d)
            else:
                fs[m] += 1  # treat as prime/undetermined
    return fs

def main():
    t0 = time.time()
    g  = 119182891324903069288022589460020572593207162963685444009526935473255725746139626528632451
    g2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    h  = 62388561396646277577754745632050284891732896280504135411098502190267395125606332503618812
    h2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    names = {'g(x9770)':g,'g2(x18274)':g2,'h(x3183)':h,'h2(x17728)':h2}
    print("=== factorizations (small primes + rho) ===", flush=True)
    for nm, v in names.items():
        fs = factor(v)
        small = {p:e for p,e in fs.items() if p < 10**9}
        print(f"  {nm} = {v}", flush=True)
        print(f"     small factors: {dict(small)}  (bits {v.bit_length()})", flush=True)
    print("\n=== pairwise gcд / relations ===", flush=True)
    import itertools
    items = list(names.items())
    for (n1,v1),(n2,v2) in itertools.combinations(items,2):
        gg = math.gcd(v1,v2)
        print(f"  gcd({n1},{n2}) = {gg}  ; {n1}/gcd={v1//gg}, {n2}/gcd={v2//gg}", flush=True)
    # h2 vs 2h, g vs g2 ratios
    print(f"\n  h2 - 2*h = {h2 - 2*h}", flush=True)
    print(f"  g - g2 = {g - g2}", flush=True)
    print(f"  gcd(g-g2, h2-2h) = {math.gcd(g-g2, h2-2*h)}", flush=True)

    # mine HUGE load constants
    A = load_atoms()
    hugeset = Counter()
    for poly in A:
        for m, c in poly.items():
            if abs(c) > 10**40: hugeset[abs(c)] += 1
    print(f"\n=== {len(hugeset)} distinct HUGE constants (|c|>1e40) ===", flush=True)
    # gcd of each huge constant with g2,h2; look for constants close to g2/h2
    hits = []
    for c in hugeset:
        for nm, q in (('g2',g2),('h2',h2),('g',g),('h',h)):
            if c % q == 0:
                hits.append((c, nm, c//q))
    print(f"HUGE constants divisible by a base quantum: {len(hits)}", flush=True)
    for c, nm, k in hits[:15]:
        print(f"   const/{nm} = {k}  (const bits {c.bit_length()})", flush=True)
    # gcd of ALL huge constants
    allg = 0
    for c in hugeset: allg = math.gcd(allg, c)
    print(f"gcd of all HUGE constants = {allg}", flush=True)
    print(f"gcd(allHuge, g2)={math.gcd(allg,g2)}, gcd(allHuge,h2)={math.gcd(allg,h2)}", flush=True)
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
