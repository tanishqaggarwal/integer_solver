"""S11 step 99: A = B = 0 HAS a solution -- eliminate w2, solve the cubic, set it.

valjac.py destroys Part XXVII's premise.  Perturbing x22152, x33462, x6418 or
x12553 moves none of w1, w2, w3, w4: those four literals are separate variables that
merely carry the same residues at this state.  The values that actually enter the
A and B are

    w1 = x12186  computed (179 free inputs move it)      w2 = x16742  FREE
    w3 = x14853  FREE                                    w4 = x24908  computed (43)
    w5 = x22162  FREE                                    w6 = x30213  FREE
    K  = x24453  the only genuine constant

so with w1, w4, w5, w6, K held, A = 0 and B = 0 are two equations in the two FREE
values w2 and w3.  Eliminate w2:

    B  =>  w = w4 - w2 = (w4+w6)(w3-w1)/(w3-w5)
    A  =>  (w5 + w1 + m + K)(m - w5)^2 = (w4+w6)^2      m = w3

a CUBIC in m.  It has exactly one root in F_p:

    w3* = 16923826268442975142014471089484050492795530131871084439458128176517372022747
    w2* = 516432665673800566800661765887332652913826924065912564591822554577222065463

and at those values A = 0 and B = 0 exactly.  **A solution to A = B = 0 exists.**
x14853 and x16742 are free variables, so the values can simply be written in; what it
costs is their own congruences a29539 (x14853 = x1308) and a26731 (x16742 = x19083),
whose targets are themselves moved by 79 and 170 free inputs -- so the cost is a
target to steer, not a wall.

Usage: close.py [state.json]
"""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from intad import jacZ
import suppfree
import fpoly as F
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'PIN_39013.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)


def report(v, tag):
    av = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
    print('%-34s score %-6d A=%s B=%s checks %s'
          % (tag, s, int(v[35389] % P == 0), int(v[6671] % P == 0), nz), flush=True)
    return s, av, nz


report(v, 'start')
w1, w2, w3, w4 = v[12186] % P, v[16742] % P, v[14853] % P, v[24908] % P
w5, w6, K = v[22162] % P, v[30213] % P, v[24453] % P
S = (w4 + w6) % P
cub = F.psub(F.pmul([(w5 + w1 + K) % P, 1],
                    [w5 * w5 % P, (-2 * w5) % P, 1]), [S * S % P])
rs = F.roots(cub)
print('\ncubic in w3 has %d root(s) in F_p' % len(rs), flush=True)
best = report(v, 'baseline')[0]
for m in rs:
    w = S * ((m - w1) % P) % P * pow((m - w5) % P, -1, P) % P
    Y1 = (w4 - w) % P
    print('\n  w3* = %d\n  w2* = %d' % (m, Y1), flush=True)
    u = list(v)
    u[14853] = (u[14853] // P) * P + m
    u[16742] = (u[16742] // P) * P + Y1
    ad.fwd(u, rounds=6)
    s, av, nz = report(u, '  after setting w3 and w2')
    print('     A = %d\n     B = %d' % (u[35389] % P, u[6671] % P), flush=True)
    # integer lift
    _, fl, SV = suppfree.build(u, modp=None)
    for _ in range(25):
        av = L.all_atom_values(u)
        todo = [a for a in range(L.NA) if a not in L.atom_out and av[a]
                and av[a] % P == 0]
        cur = L.NEQ - len(L.failing_eqs(av))
        moved = False
        for c in todo:
            mm = suppfree.atom_supp(c, u, SV, modp=None)
            for i in range(len(fl)):
                if not ((mm >> i) & 1):
                    continue
                t = fl[i]
                g = jacZ(t, u, [c]).get(c, 0)
                if not g or g % P or av[c] % g:
                    continue
                w2 = list(u)
                w2[t] = w2[t] - av[c] // g
                ad.fwd(w2, rounds=6)
                a2 = L.all_atom_values(w2)
                if a2[c] == 0 and L.NEQ - len(L.failing_eqs(a2)) >= cur:
                    u, av, moved = w2, a2, True
                    break
            if moved:
                break
        if not moved:
            break
    s, av, nz = report(u, '  after the integer lift')
    T.save(u, os.path.join(HERE, 'CLOSE_%d.json' % s))
    print('  saved CLOSE_%d.json' % s, flush=True)
    if s > best:
        best = s
print('\nbest %d' % best)
