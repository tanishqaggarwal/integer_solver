import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import atomlib as A
p = A.p
BASE = '/home/user/integer_solver/solve_lab'
v = A.load_json(BASE + '/best_agentA_39022.json')

for leaf in [4432, 7068, 17325, 9413, 19964, 2099]:
    print(f"\n===== x_{leaf} appears in atoms: =====")
    for ai in A.VAR_ATOMS[leaf]:
        val = A.eval_atom(ai, v)
        eqs = A.ATOM_EQS[ai]
        print(f"  atom {ai}: {A.ATOM_REPR[ai]!r}")
        print(f"      val%p={'0' if val%p==0 else 'NONZERO'}  vars={sorted(A.ATOM_VARS[ai])}  in {len(eqs)} eqs: {eqs[:12]}")
