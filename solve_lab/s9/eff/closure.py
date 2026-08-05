"""For each congruence-breaking variable, take the FULL closure of its block
(iterate S <- S u eqs(confined atoms + the breaker's extra atoms)) and price it."""
import pickle, sys, time
import lib as L, model as MD, opt

v0 = opt.init()
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])


def closure(S, seed_atoms, rounds=6, cap=400):
    S = set(S)
    A = set(seed_atoms)
    for _ in range(rounds):
        for a in list(A):
            S |= set(L.atom2eq.get(a, ()))
        if len(S) > cap:
            return None
        newA = set(MD.confined_atoms(S))
        # also pull in every atom that occurs in S and is nearly confined
        if newA <= A and all(set(L.atom2eq.get(a, ())) <= S for a in A):
            A = newA
            break
        A |= newA
    return frozenset(S)


if __name__ == '__main__':
    rows = pickle.load(open('breakers.pkl', 'rb'))
    A13 = set(MD.confined_atoms(S13))
    print('breaker closures:')
    for cost, x, b1, b2, ea, ce in rows[:20]:
        S = closure(S13, A13 | set(ea))
        if S is None:
            print(f'  x_{x}: closure exceeds cap')
            continue
        mod = MD.build(S, v0, verbose=False)
        D = opt.rank_of(mod)
        print(f'  x_{x:<6d} closed |S|={len(S):<4d} |A|={len(mod["A"]):<4d} knobs={len(mod["knobs"]):<4d} '
              f'D={D:<4d} f={len(S)-D:<4d} minfail>={len(S)-D+1}')
