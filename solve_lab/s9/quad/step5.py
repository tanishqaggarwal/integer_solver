"""Step 5: check the couplings of x_16742 / x_12186 / x_6418 / x_12553 / x_22162 / x_30213."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

v = H.load_assignment('quad/stateA1.json')
for u in (16742, 19083, 9254, 12186, 30454, 26865, 33612, 6418, 12553, 22162, 30213, 7497,
          11436, 17286, 13632, 8386, 21868, 22820, 14393):
    a = definer.get(u)
    print(f'x_{u:<6d} def={a}  free={a is None}  val_bits={v[u].bit_length()}  occurs_in={len(var_atoms[u])} atoms')
    if a is not None:
        print(f'          {src[a][:140]}')

print('\n# which CHECK atoms mention each candidate var')
for u in (16742, 12186, 6418, 12553, 22162, 30213, 30454, 26865):
    ca = [a for a in var_atoms[u] if a in checkset]
    ga = [a for a in var_atoms[u] if a not in checkset]
    print(f'x_{u}: checks={ca}  gates={len(ga)}')
    for a in ca:
        print(f'    [{a}] resid={"NZ" if evalpoly(resid_poly[a],v)!=0 else "0 "}  {src[a][:130]}')
