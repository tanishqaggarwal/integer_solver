"""S11 step 74: CONSTRUCTIVE gadget closing -- stop solving, start substituting.

Every failing check in the canonical frame has the same shape:

    a21617 = 11436039*(x14623 - x27522) - p*x5040        x5040  free (handle)
    a29539 = 12846437*(x14853 - x1308 ) - p*x30163       x30163 free (handle)

x14623 and x14853 are FREE, and the other side's free-input support (11 and 79
inputs) contains neither of them -- so the two sides are independent.  That makes
the fix a substitution, not a search:

    x14623 <- x14623 - ((x14623 - x27522) mod p)   =>  p | c*(x14623 - x27522)
    x5040  <- c*(x14623 - x27522)/p                =>  the atom is exactly zero

both steps exact over Z.  All the linear-algebra vetoes in this lab priced MOVES;
none of them priced this, because it is not a move in the tangent space -- it is a
jump of a full residue class.

Usage: gfix.py [state.json] [out.json]
"""
import os, sys, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'mod9118_0.json'
out = sys.argv[2] if len(sys.argv) > 2 else 'gfix.json'
FREE = set(t for t in range(L.NVARS) if t not in L.definer)

v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def report(tag):
    av = L.all_atom_values(v)
    bad = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    gb = [a for a in L.atom_out if av[a]]
    s = L.NEQ - len(L.failing_eqs(av))
    print('%-34s score %d   failing checks %s   broken gates %d'
          % (tag, s, bad, len(gb)), flush=True)
    return s, bad


report('start')
# (coef, freeside, otherside, handle)
GAD = [(11436039, 14623, 27522, 5040),
       (12846437, 14853, 1308, 30163)]
for c, lhs, rhs, h in GAD:
    print('\n--- gadget  %d*(x%d - x%d) - p*x%d ---' % (c, lhs, rhs, h))
    print('    x%d free: %s;  handle x%d free: %s'
          % (lhs, lhs in FREE, h, h in FREE))
    d = (v[lhs] - v[rhs]) % P
    print('    (x%d - x%d) mod p = %d' % (lhs, rhs, d))
    v[lhs] -= d
    ad.fwd(v, rounds=6)
    num = c * (v[lhs] - v[rhs])
    print('    after the jump: p | c*(A-B)? %s' % (num % P == 0))
    if num % P:
        print('    *** the other side MOVED with the jump -- not independent')
        continue
    v[h] = num // P
    ad.fwd(v, rounds=6)
    report('after closing this gadget')

s, bad = report('\nFINAL')
T.save(v, os.path.join(HERE, out))
print('saved %s' % out)
