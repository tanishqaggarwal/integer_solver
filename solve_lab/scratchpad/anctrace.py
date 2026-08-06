import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

# find the c-side pin: an atom of form coef*(x_24908 - x_?) - slack, and x_12186 pin (x_23927)
print("=== atoms containing x_24908 (looking for its pin) ===")
for ai in A.VAR_ATOMS[24908]:
    if len(A.ATOM_VARS[ai]) <= 6:
        print(f"  atom {ai}: {A.ATOM_REPR[ai][:90]!r} vars={sorted(A.ATOM_VARS[ai])}")
print("=== atoms containing x_12186 (small) ===")
for ai in A.VAR_ATOMS[12186]:
    if len(A.ATOM_VARS[ai]) <= 6:
        print(f"  atom {ai}: {A.ATOM_REPR[ai][:90]!r} vars={sorted(A.ATOM_VARS[ai])}")

# free-input ancestors
key = {'x_1308':1308,'x_23927':23927,'x_19083':19083,'x_14853':14853,
       'x_12186':12186,'x_16742':16742,'x_24908':24908}
anc = {}
for name,x in key.items():
    a = H.anc.get(x, set())
    anc[name] = a
    print(f"\n{name} (x_{x}): {'FREE' if x in H.freeinp else 'GATE'}, #free-ancestors={len(a)}")
    if len(a) <= 20:
        print(f"   ancestors: {sorted(a)}")

# overlaps for the key consistency conditions
print("\n=== overlap x_1308 vs x_23927 (need x_1308 == x_23927 mod p) ===")
print(f"  |anc(1308)|={len(anc['x_1308'])}, |anc(23927)|={len(anc['x_23927'])}, shared={len(anc['x_1308'] & anc['x_23927'])}")
print(f"  1308-only: {sorted(anc['x_1308'] - anc['x_23927'])[:20]}")
print(f"  23927-only: {sorted(anc['x_23927'] - anc['x_1308'])[:20]}")
