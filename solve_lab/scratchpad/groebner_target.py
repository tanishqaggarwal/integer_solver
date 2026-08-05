import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

RES = [6221,18868,32628,41908,42325,43165,43419,43447,44713,45258]
allv = set()
for ai in RES: allv |= A.ATOM_VARS[ai]
print(f"vars across 10 atoms: {len(allv)}")
free = sorted(v for v in allv if v in H.freeinp)
gates = sorted(v for v in allv if v not in H.freeinp)
print(f"  free inputs directly: {len(free)}: {free}")
print(f"  gate outputs: {len(gates)}")

# transitive free-input ancestors of all these atoms' gate outputs
anc_all = set(free)
for g in gates:
    anc_all |= H.anc.get(g, set())
print(f"\ntransitive free-input ancestors (true unknowns): {len(anc_all)}")

# which are bits?
import json
bits = set(json.load(open('all_bits.json')))
anc_bits = anc_all & bits
anc_nonbit = anc_all - bits
print(f"  bits (0/1): {len(anc_bits)}")
print(f"  non-bit free inputs: {len(anc_nonbit)}: {sorted(anc_nonbit)}")

# Evaluate the 10 atoms at fullcore_fix and at 39022
for name, path in [('fullcore_fix','fullcore_fix.json'), ('39022','best_agentA_39022.json')]:
    v = A.load_json(path)
    nz = [(ai, A.eval_atom(ai, v) % p) for ai in RES]
    nzc = sum(1 for ai, val in nz if val != 0)
    print(f"\n{name}: {nzc}/10 nonzero mod p: {[(ai) for ai,val in nz if val!=0]}")
