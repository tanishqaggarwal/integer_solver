import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
B.regime11()
F=set(H.fails())
print('current fails:',sorted(F))
# knobs and the eqs they appear in (via eqvars)
knobs=[8731,9118,950,6947,33168,9413,17325,4432,7068]
# For a free input, which eqs contain it? eqvars are per-equation var sets.
for kn in knobs:
    eqs=[i for i,vs in enumerate(H.eqvars) if kn in vs]
    outside=[e for e in eqs if e not in F]
    print(f'x_{kn}: in {len(eqs)} eqs; {len(outside)} OUTSIDE current fails: {outside[:15]}')
