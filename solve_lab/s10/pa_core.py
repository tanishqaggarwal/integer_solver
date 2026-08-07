"""Forced-core bound.

Ca = x_7068 - x_2099 is nonzero mod P, and a22229 is the ONLY atom containing both
x_7068 and x_2099.  a22229 = 0 forces x_7068-x_2099 = 7376877*x_642, and a35762 = 0
forces x_642 = P*x_17325, hence Ca == 0 (mod 7376877*P).  So any placement that does not
move Ca must contain a22229 or a35762 -- both with a 9-equation footprint.

For each forced core we exhaustively search all extensions by <=3 equations for the best
deficiency |E| - |A(E)| (primitive atoms only)."""
import os, sys, itertools, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
PRIM=[a for a in range(L.NA) if len(EQS[a])>1]
both=[a for a in range(L.NA) if 7068 in L.avars[a] and 2099 in L.avars[a]]
print('atoms containing both x_7068 and x_2099:',both)
print('atoms containing x_642:',sorted(L.var_atoms[642]))
for core in (22229,35762,22230):
    E=set(EQS[core])
    miss={b:(EQS[b]-E) for b in PRIM if len(EQS[b]-E)<=4}
    pool=sorted(set().union(*[m for m in miss.values() if 1<=len(m)<=3]))
    base=[b for b,m in miss.items() if not m]
    print(f'\ncore a{core}: |E|={len(E)} A(E)={base} deficiency={len(E)-len(base)}  pool={len(pool)}')
    for k in (1,2,3):
        best=None
        for X in itertools.combinations(pool,k):
            Xs=set(X)
            n=sum(1 for b,m in miss.items() if m<=Xs)
            d=len(E)+k-n
            if best is None or d<best[0]: best=(d,X,n)
        print(f'   best extension |X|={k}: deficiency={best[0]} (|E|={len(E)+k}, |A|={best[2]}) X={best[1]}')
