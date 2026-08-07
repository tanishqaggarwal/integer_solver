"""S11 step 72: decompile the CANONICAL-frame residual.

Step 71 showed the residual atoms are one multiplication each, and that almost every
operand chain collapses onto x26064 = p.  The canonical frame is where forward
evaluation makes every gate hold by construction, so there the residual is only the
CHECK atoms -- a much smaller object than the seven of the witness frame.  Print the
failing checks of a given state as polynomials, then unfold every operand until it
is either free or a constant, so the residual can be read as an arithmetic identity
instead of attacked as a linear system.

Usage: decomp2.py [state.json] [DEPTH]
"""
import os, sys, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'mod9118_0.json')
DEPTH = int(sys.argv[2]) if len(sys.argv) > 2 else 40
v = L.load(path if os.path.isabs(path) else os.path.join(HERE, path))
FREE = set(t for t in range(L.NVARS) if t not in L.definer)


def poly(a):
    ts = L.polys[a]
    it = ts.items() if isinstance(ts, dict) else ts
    out = []
    for mono, c in it:
        if isinstance(mono, tuple) and mono:
            out.append('%+d*%s' % (c, '*'.join('x%d' % m for m in mono)))
        else:
            out.append('%+d' % c)
    return ' '.join(out)


av = L.all_atom_values(v)
CH = [a for a in range(L.NA) if a not in L.atom_out]
BAD = [a for a in CH if av[a]]
GBAD = [a for a in L.atom_out if av[a]]
print('state %s: score %d' % (os.path.basename(path), L.NEQ - len(L.failing_eqs(av))))
print('  failing CHECK atoms: %s' % BAD)
print('  broken GATE atoms:   %d %s' % (len(GBAD), GBAD[:12]))
print()
for a in BAD:
    print('=== a%d = %s   (value %d bits, = %s mod p)'
          % (a, poly(a), av[a].bit_length(), av[a] % P))
    seen, front = set(), sorted(L.avars[a])
    for lv in range(DEPTH):
        nxt = set()
        for t in sorted(front):
            if t in seen:
                continue
            seen.add(t)
            if t in FREE:
                print('    x%-6d FREE  (%d bits)%s'
                      % (t, v[t].bit_length(),
                         '  == 0 mod p' if v[t] % P == 0 else ''))
                continue
            d = L.definer[t]
            print('    x%-6d <- a%-6d %s' % (t, d, poly(d)[:130]))
            for w in L.avars[d]:
                if w != t and w not in seen:
                    nxt.add(w)
        front = sorted(nxt)
        if not front:
            break
    print('    (%d vars, %d free)\n' % (len(seen), sum(1 for t in seen if t in FREE)))
