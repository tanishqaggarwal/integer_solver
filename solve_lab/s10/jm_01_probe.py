"""jm step 1: exact broken-CHECK sets for the key single moves, and the
special values of the branch inputs x_2081 / x_4287."""
import os, sys, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
w = J.base_state()
avb = L.all_atom_values(w)
R0 = J.resid(w)
print('base out12/score', J.cost(w)[:2], 'resid moved-from', [str(x)[:12] for x in R0])


def probe(moves, tag='', show=True):
    v = list(w)
    for u, val in moves.items():
        v[u] = val
    J.fwd2(v, 2)
    c, s, f, av = J.cost(v)
    nz = [a for a in range(L.NA) if av[a] and a not in J.SS]
    R = J.resid(v)
    mv = ''.join(n for n, a, b in zip('ABC', R, R0) if a != b)
    if show:
        print(f'  {tag:<44} out12={c:<4} score={s:<6} moves={mv:<3} '
              f'brokenchecks={nz[:14]}', flush=True)
    return c, s, nz, mv, v


if __name__ == '__main__':
    print('\n--- single movers, exact broken-check sets ---')
    for u in (6418, 24548, 12553, 14853, 4287, 2081, 13195, 28730, 7068):
        probe({u: w[u] + 1000003}, f'x_{u} += 1000003')
    print('\n--- branch inputs at their special values ---')
    for u, val in ((2081, 0), (2081, 2), (2081, -1), (4287, 0), (4287, 1), (4287, -1)):
        probe({u: val}, f'x_{u} = {val}')
    print('\n--- x_2081=0 then move x_6418 (a3576 should go dead) ---')
    for d in (1, 1000003, P):
        probe({2081: 0, 6418: w[6418] + d}, f'x_2081=0, x_6418 += {str(d)[:12]}')
    print('\n--- pure p-multiple moves (no residue change) ---')
    for u in (6418, 24548, 12553, 14853):
        probe({u: w[u] + P}, f'x_{u} += p')
