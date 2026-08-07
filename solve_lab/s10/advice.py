"""S11 step 85: THE FREE CONTENT OF THE INSTANCE IS THIRTEEN 296-BIT NUMBERS.

Of the 7,273 free inputs, 7,252 are zero at every good state in this lab, a handful
are the large values build7.py writes, and exactly THIRTEEN carry 296 bits:

    x6418  x8778  x12553 x14623 x14853 x16742 x22152 x22162 x22649 x24548
    x30213 x31339 x33462

That is the whole search space.  296 bits is 256 + 40: a field element plus ~40 bits
of k*p slack, which is what a gadget `c*(x - B) - p*h` wants -- x may be any member
of B's residue class, and h absorbs the difference.

So the instance is thirteen advice values, eleven of which are already right; the two
that are wrong are exactly the two remaining congruences

    x24548 ≡ x25442 (mod p)        x14853 ≡ x1308 (mod p)

This module prints, for each of the thirteen, every atom it occurs in, which of them
are checks, what congruence each check imposes on it, and whether that congruence
currently holds -- so the constraints on the advice can be read instead of searched.

Usage: advice.py [state.json]
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'B7_39004.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
av = L.all_atom_values(v)
FREE = [t for t in range(L.NVARS) if t not in L.definer]
ADV = [t for t in FREE if 280 <= v[t].bit_length() <= 330]
print('%s: score %d;  %d advice inputs'
      % (src, L.NEQ - len(L.failing_eqs(av)), len(ADV)), flush=True)


def poly(a):
    ts = L.polys[a]
    it = ts.items() if isinstance(ts, dict) else ts
    return ' '.join(('%+d*%s' % (c, '*'.join('x%d' % m for m in mono)))
                    if isinstance(mono, tuple) and mono else '%+d' % c
                    for mono, c in it)


for t in ADV:
    q, r = divmod(v[t], P)
    print('\n=== x%-6d = %d*p + r   (r has %d bits, k has %d bits)'
          % (t, q, r.bit_length(), q.bit_length()))
    for a in sorted(L.var_atoms[t]):
        kind = 'GATE ->x%d' % L.atom_out[a][1] if a in L.atom_out else 'CHECK'
        flag = '' if av[a] == 0 else '   <<< NONZERO'
        print('    a%-6d %-12s %s%s' % (a, kind, poly(a)[:120], flag))
print('\n--- what each advice input is congruent to (from its two-sided gadgets) ---')
for t in ADV:
    for a in sorted(L.var_atoms[t]):
        if a in L.atom_out:
            continue
        pl = L.polys[a]
        if len(pl) > 6:
            continue
        terms = [(m, c) for m, c in pl.items()]
        lin = [(m[0], c) for m, c in terms if len(m) == 1]
        if len(lin) == 2 and any(w == t for w, _ in lin):
            (w1, c1), (w2, c2) = lin
            other = w2 if w1 == t else w1
            co = c1 if w1 == t else c2
            d = (v[t] - v[other]) % P
            print('    x%-6d ≡ x%-6d (mod p) via a%-6d coefficient %d : %s'
                  % (t, other, a, co, 'HOLDS' if d == 0 else 'FAILS (diff %d)' % d))
