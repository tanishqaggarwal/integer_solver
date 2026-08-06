#!/usr/bin/env python3
"""Exploit the kernel algebraically. Around the x_3368 slack:
  a27973: x_9770 = x_35186 + x_3368
  a1660 : x_3368 = x_12779 * x_24026
  a1657 : x_12779 + x_14402 = 1        => x_14402 = 1 - x_12779
  a1813 : x_14402 * x_24026 = 321447 * x_38215
Adding a1660 and a1813 (t=x_12779):  t*u + (1-t)*u = u = (x_9770-x_35186) + 321447*x_38215.
If x_38215 is a fixed constant, u=x_24026 and t=x_12779 are DETERMINED by x_9770 (=x_18274 at witness),
x_35186, x_38215.  Analogous chain for x_3183/x_10466/x_27116.  Compute all constants & check whether
the witness's t must be a small integer achievable by x_12779=x_23380*x_36336.
"""
import json
from math import gcd
from confluent_eval5 import build5, make_forward
from propagate import load_atoms, atom_vars, NVARS
from collections import deque

def main():
    atoms = load_atoms()
    best = json.load(open('best/best_partial_39019.json')); bv = {int(k[2:]): v for k, v in best.items()}
    def val(v): return bv.get(v, 0)
    control = set(json.load(open('control_bits.json')))

    # confirm the atoms
    def showatom(a):
        return {m: c for m, c in atoms[a].items()}
    print("a27973:", showatom(27973))
    print("a1660 :", showatom(1660))
    print("a1657 :", showatom(1657))
    print("a1813 :", showatom(1813))
    print("a1815 :", showatom(1815))
    print("a27978:", showatom(27978))
    print("a27976:", showatom(27976))

    # is x_38215 constant? build its cone (via all atoms, treating each var's shortest atom as def)
    # simpler: does x_38215 depend on any control bit? use forward-eval single-flip
    Aat, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    base = solve(list(bestval), [])
    for w in (38215, 29437, 35186, 1642, 4028, 27787, 32103):
        movers = [b for b in control if solve(list(bestval), [b])[w] != base[w]]
        print(f"x_{w} = {base[w]}  (moved by {len(movers)} single bits)")

    g = 119182891324903069288022589460020572593207162963685444009526935473255725746139626528632451
    D1 = base[9770] - base[18274]
    print(f"\nD1 = x_9770-x_18274 (best) = {D1}")
    print(f"x_35186(best) = {base[35186]}  (= m*g? m={base[35186]//g if base[35186]%g==0 else 'NO'})")
    print(f"x_38215(best) = {base[38215]}")
    print(f"321447*x_38215 = {321447*base[38215]}")
    # bridge: need x_3368 = x_18274 - x_35186 at the (merged) solution. At best x_9770=x_35186 (x_3368=0).
    # For witness x_9770=x_18274: x_3368 = x_18274 - x_35186.
    # u = x_24026 = (x_18274 - x_35186) + 321447*x_38215 ; t=x_12779 = (x_18274-x_35186)/u
    # But x_18274, x_35186 depend on bits. Explore: if x_38215 const C and x_35186=m*g:
    C = base[38215]
    print(f"\nIf x_38215 == {C} (const) and x_35186=m*g, then for target x_18274=V (=x_9770):")
    print(f"  x_24026 = (V - m*g) + 321447*{C}")
    print(f"  x_12779 = (V - m*g) / x_24026  must be a small integer (achievable x_23380*x_36336).")
    # gcd(g, 321447)
    print(f"gcd(g,321447)={gcd(g,321447)}; 321447 factors:", [p for p in range(2,321448) if 321447%p==0][:6])

if __name__ == '__main__':
    main()
