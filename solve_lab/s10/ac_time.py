"""S12 step 2: timing probe for the expensive primitives."""
import os, sys, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
FREE = set(ad.FREE)
random.seed(1)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
t=time.time(); av0 = L.all_atom_values(v0); print(f'all_atom_values {time.time()-t:.2f}s')
t=time.time(); f = L.failing_eqs(av0); print(f'failing_eqs {time.time()-t:.2f}s -> {len(f)}')
vm0=[x%P for x in v0]
t=time.time(); g=ad.grad(21617, vm0); print(f'grad {time.time()-t:.2f}s supp {len(g)}')
CHECKS=sorted(a for a in range(L.NA) if a not in L.atom_out)
print('nchecks', len(CHECKS))
t=time.time(); c=jac_column(list(g)[0], v0, vm0, CHECKS); print(f'jac_column {time.time()-t:.2f}s supp {len(c)}')
z = 12054
v=list(v0); v[z]=random.randrange(1,1<<48)
t=time.time(); ad.fwd(v, rounds=6); print(f'ad.fwd r6 {time.time()-t:.2f}s')
ch = {i for i in range(L.NVARS) if v[i]!=v0[i]}
print('changed vars', len(ch))
t=time.time(); ta = L.touched_atoms(v, av0, ch); print(f'touched_atoms {time.time()-t:.2f}s -> {len(ta)}')
v2=list(v0); v2[z]=v[z]
t=time.time(); ch2,_ = L.ripple(v2, {z: v[z]}); print(f'L.ripple {time.time()-t:.2f}s changed {len(ch2)}')
diff = sum(1 for i in range(L.NVARS) if v[i]!=v2[i])
print('fwd vs ripple differing vars:', diff)
t=time.time(); ad.fwd(v2, rounds=2); print(f'ripple+fwd2 {time.time()-t:.2f}s')
print('after: differing', sum(1 for i in range(L.NVARS) if v[i]!=v2[i]))
