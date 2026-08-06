import sys, os, json
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p
RES = [6221,18868,32628,41908,42325,43165,43419,43447,44713,45258]

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v] = vA.get(v, 0)
H.forward()
base = {ai: A.eval_atom(ai, H.val) % p for ai in RES}

nonbit = [124, 371, 893, 2263, 3067, 3192, 3203, 3272, 3325, 3473, 4210, 4715, 5616, 6083, 6882, 7256, 7583, 7994, 8055, 8060, 8183, 8512, 8605, 8971, 9776, 10261, 11005, 11049, 11080, 11559, 11739, 12104, 12398, 13638, 13771, 14375, 14412, 15123, 15692, 16259, 16441, 16773, 16787, 17225, 18991, 19093, 19275, 19450, 19714, 19844, 20007, 20555, 20824, 21815, 21846, 22321, 22727, 23125, 23971, 25012, 26040, 26348, 26489, 26602, 26813, 27156, 27305, 27616, 28246, 28486, 28548, 29261, 30060, 30176, 30454, 30468, 31339, 31460, 31729, 32203, 32881, 33129, 33229, 33708, 36804, 36925, 37589, 38101, 38258, 38460, 38480]

couplers = {}
for m in nonbit:
    old = H.val[m]
    H.val[m] = old + 1
    H.forward()
    moved = [ai for ai in RES if A.eval_atom(ai, H.val) % p != base[ai]]
    if moved:
        couplers[m] = moved
    H.val[m] = old
H.forward()

print(f"non-bit inputs that move >=1 of the 10 atoms: {len(couplers)}")
for m in sorted(couplers, key=lambda x: -len(couplers[x])):
    print(f"  x_{m}: moves atoms {couplers[m]}")

# also test the 256 bits quickly
bits = sorted(set(json.load(open('all_bits.json'))))
bit_couplers = {}
for m in bits:
    old = H.val[m]
    H.val[m] = 1 - old if old in (0,1) else old + 1
    H.forward()
    moved = [ai for ai in RES if A.eval_atom(ai, H.val) % p != base[ai]]
    if moved: bit_couplers[m] = moved
    H.val[m] = old
H.forward()
print(f"\nBITS that move >=1 of the 10 atoms: {len(bit_couplers)}")
print(f"  sample: {dict(list(bit_couplers.items())[:10])}")
