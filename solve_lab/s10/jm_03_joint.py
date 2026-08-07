"""jm step 3: the TRACKING joint moves.

a29539 = 12846437*(x_14853 - x_1308) - x_29967   (x_29967 = p * x_30163)
a7930  =  9367949*(x_24548 - x_25442) - x_7927   (x_7927  = p * x_11052)
a21617 = 11436039*(x_14623 - x_27522) - x_36864  (x_36864 = p * x_5040)
a3576  = x_2081*(x_6418 - C) - 15804267*x_26777  (x_26777 = p * x_3387)

Moving a driver breaks the pin atom; move the FREE side of the pin to track it
exactly and the pin closes again.  Price every such tracking move.
"""
import os, sys, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
w = J.base_state()
R0 = J.resid(w)
avb = L.all_atom_values(w)

# (pin atom, free tracker var, tracked var, integer handle var, handle free input)
PINS = [(29539, 14853, 1308, 29967, 30163),
        (7930, 24548, 25442, 7927, 11052),
        (21617, 14623, 27522, 36864, 5040),
        (3576, None, None, 26777, 3387)]


def report(v, tag):
    c, s, f, av = J.cost(v)
    nz = [a for a in range(L.NA) if av[a] and a not in J.SS]
    R = J.resid(v)
    mv = ''.join(n for n, a, b in zip('ABC', R, R0) if a != b)
    print(f'  {tag:<52} out12={c:<4} score={s:<6} moves={mv:<4} broken={nz}',
          flush=True)
    return c, s, nz, mv


def track(v, pin, tracker, tracked):
    """set the free tracker so the pin closes again, after v was moved."""
    v[tracker] = v[tracker] + (v[tracked] - w[tracked])
    J.fwd2(v, 2)
    return v


if __name__ == '__main__':
    print('handles (should be p-quantised):')
    for a, tr, td, h, hf in PINS:
        print(f'  a{a}: handle x_{h} = {L.atom_src[J.definer[h]][:70]}')
    print()

    print('--- J1: move x_7068, track x_1308 with x_14853  (congruence 1) ---')
    for d in (1, 3, 1000003):
        v = list(w); v[7068] += d; J.fwd2(v, 2)
        dd = v[1308] - w[1308]
        v[14853] += dd; J.fwd2(v, 2)
        report(v, f'x_7068+={d}, x_14853+={str(dd)[:18]}')

    print('--- J1b: move x_6418 instead (congruence 1 by the other side) ---')
    for d in (1, 1000003):
        v = list(w); v[6418] += d; J.fwd2(v, 2)
        report(v, f'x_6418+={d}')

    print('--- J2: move x_28730, track x_25442 with x_24548 (congruence 2) ---')
    for d in (1, 3, 1000003):
        v = list(w); v[28730] += d; J.fwd2(v, 2)
        dd = v[25442] - w[25442]
        v[24548] += dd; J.fwd2(v, 2)
        report(v, f'x_28730+={d}, x_24548+={str(dd)[:18]}')

    print('--- J2b: + also track x_27522 with x_14623 (kills a21617) ---')
    for d in (1, 1000003):
        v = list(w); v[28730] += d; J.fwd2(v, 2)
        v[24548] += v[25442] - w[25442]
        v[14623] += v[27522] - w[27522]
        J.fwd2(v, 2)
        report(v, f'x_28730+={d} + 2 trackers')

    print('--- J3: BOTH congruences at once ---')
    for d1, d2 in ((1, 1), (1000003, 1000033)):
        v = list(w); v[7068] += d1; v[28730] += d2; J.fwd2(v, 2)
        v[14853] += v[1308] - w[1308]
        v[24548] += v[25442] - w[25442]
        v[14623] += v[27522] - w[27522]
        J.fwd2(v, 2)
        report(v, f'x_7068+={d1}, x_28730+={d2}, 3 trackers')
