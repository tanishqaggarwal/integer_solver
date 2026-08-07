"""S12 step 17: after the repair engine wins the collateral back, is the
ACTIVATION still live?  If the repaired state scores 39009 but has a larger
cluster gradient support, the knobs were bought for free."""
import os, sys, glob, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad, ac_lib as A
P = ad.P
B = A.Base(os.path.join(HERE,'mod9118_0.json'))
BAD = [21617, 29539]
supp0 = A.grad_supp(B.v0, BAD)
print(f'base: score {B.score0}  nonzero {sorted(B.nz0)}  cluster support {len(supp0)}')
for f in sorted(glob.glob(os.path.join(HERE,'engine_ac_*.json'))):
    v = L.load(f)
    av = L.all_atom_values(v)
    nz = sorted(a for a in range(L.NA) if av[a])
    sc = L.NEQ - len(L.failing_eqs(av))
    nzchk = [a for a in nz if a in A.CHECKSET]
    supp = A.grad_supp(v, BAD)
    diff = sum(1 for i in range(L.NVARS) if v[i] != B.v0[i])
    nzf = sum(1 for u in ad.FREE if v[u])
    print(f'{os.path.basename(f)}: score {sc}  nonzero {nz}\n'
          f'    cluster support {len(supp)} (base {len(supp0)}, new knobs '
          f'{len(supp-supp0)}, lost {len(supp0-supp)})  vars differing from base {diff}  '
          f'nonzero free inputs {nzf}')
