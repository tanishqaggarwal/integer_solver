"""S11 step 73: read each failing CHECK as a modular-reduction gadget.

In the canonical frame every gate holds by construction, so the entire problem is
"choose the free inputs so that every CHECK atom vanishes".  At mod9118_0 only four
checks are nonzero -- four constraints against 7,273 free inputs.  And a21617 reads

    11436039*x14623 - 11436039*x27522 - x36864 ,   x36864 = x986*x5040,  x986 = p

i.e. c*(A - B) - p*handle:  it asserts A == B (mod p) with x5040 absorbing the
quotient.  If the handle is free and used nowhere else, the check is satisfiable the
instant p divides c*(A-B); if A or B has a free summand, we can even force that.
So the real question for each failing check is not "can a linear solver move it"
but "which of its operands are free, and what else do they touch".

Usage: gadget.py [state.json]
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
path = sys.argv[1] if len(sys.argv) > 1 else 'mod9118_0.json'
v = L.load(path if os.path.isabs(path) else os.path.join(HERE, path))
FREE = set(t for t in range(L.NVARS) if t not in L.definer)
av = L.all_atom_values(v)


def poly(a):
    ts = L.polys[a]
    it = ts.items() if isinstance(ts, dict) else ts
    return ' '.join(
        ('%+d*%s' % (c, '*'.join('x%d' % m for m in mono)))
        if isinstance(mono, tuple) and mono else '%+d' % c for mono, c in it)


def usage(t):
    """Everything x_t touches: atoms it appears in, and how they are used."""
    ats = sorted(L.var_atoms[t])
    gates = [a for a in ats if a in L.atom_out]
    checks = [a for a in ats if a not in L.atom_out]
    eqs = set()
    for a in ats:
        eqs |= set(L.atom2eq[a])
    return ats, gates, checks, sorted(eqs)


BAD = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
print('state %s: score %d; failing checks %s'
      % (os.path.basename(path), L.NEQ - len(L.failing_eqs(av)), BAD))
print()
for a in BAD:
    print('=== a%d = %s' % (a, poly(a)))
    print('    value mod p = %d' % (av[a] % P))
    for t in sorted(L.avars[a]):
        ats, gates, checks, eqs = usage(t)
        tag = 'FREE' if t in FREE else 'a%d' % L.definer[t]
        print('    x%-6d %-8s %4d bits  |  appears in %d atoms (%d gate, %d check),'
              ' %d equations%s'
              % (t, tag, v[t].bit_length(), len(ats), len(gates), len(checks),
                 len(eqs), '   <== p' if v[t] == P else ''))
        if t in FREE and len(ats) <= 4:
            for b in ats:
                print('           a%-6d %s%s' % (b, poly(b)[:110],
                                                 '' if b in L.atom_out else '  [CHECK]'))
    print()
